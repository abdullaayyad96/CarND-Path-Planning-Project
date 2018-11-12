#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <limits>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

using namespace std;

typedef std::numeric_limits< double > dbl;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/2)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;
	
	while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1)))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);
	
	//cout << s << "\t" << d << "\t" << x << "\t" << y << endl;

	return {x,y};

}

vector<vector<double>> potential_funtion(double start_s, double start_d, double start_x, double start_y, double start_yaw, double cur_speed, double ref_velocity_s, double max_acceleration, int start_point, vector<vector<double>> sensor_fusion_sd_frame, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	vector<vector<double>> xy_points;
	//auto start_xy = getXY(start_s, start_d, maps_s, maps_x, maps_y);
	//cout << "start" << "\t" << start_s << "\t" << start_d << "\t" << start_xy[0] << "\t" << start_xy[1] << endl;
	double s = start_s;
	double d = start_d;
	double x = start_x;
	double y = start_y;
	double speed_x = cur_speed * cos(start_yaw);
	double speed_y = cur_speed * sin(start_yaw);

	//Construct potential function
	//PF =  - ref_velocity * s + w_b1 * exp(- w_b2 * d^2) +  w_b1 * exp(- w_b2 * (d-12)^2) +  summation of w_v1 exp(- w_v2 * (s - s_v)^2 - w_v3 *(d - d_v)^2 ) for every vehicle
	double w_b1 = 2000;
	double w_b2 = 1;
	double w_v1 = 100;
	double w_v2 = 0.01;
	double w_v3 = 10;
	double w_v4 = 0.2;
	double w_lane = 0;

	double speed_s = speed_x;

	double new_x, new_y, new_yaw, new_speed_x, new_speed_y, acc_x, acc_y, speed_mag, acc_mag;
	vector<double> new_points, new_sd;

	for (int i = start_point; i < 50; i++)
	{
		//iterate through data points
		
		//calculate gradient of PF
		double PF_grad_s = -cur_speed;
		double PF_grad_d = 0;

		if (cur_speed < ref_velocity_s)
			cur_speed += 0.25;
		
		PF_grad_d += -2 * w_b1 * w_b2 * d * exp(-w_b2 * d*d) - 2 * w_b1 * w_b2 * (d - 12) * exp(-w_b2 * pow(d - 12, 2));

		//quadratic function to keep lane
		int lane = (int)d / 4;
		double d_lane = 2 + 4 * lane;
		PF_grad_d += 2 * w_lane * (d - d_lane);
		

		
		for (int j = 0; j < sensor_fusion_sd_frame.size(); j++)
		{
			//iterate through  cars

			//linearly interpolate the position of each car
			double s_v = sensor_fusion_sd_frame[j][0] + i * sensor_fusion_sd_frame[j][2] * 0.02;
			double d_v = sensor_fusion_sd_frame[j][1];// +i * sensor_fusion_sd_frame[j][3] * 0.02;

			while (s_v > 6945.554)
				s_v -= 6945.554;

				//update PF
				//if(abs(d-d_v)<2)
			if ( ((s - s_v) < 5)  && ( (s - s_v) > -20 ) ) {
				PF_grad_s += 2 * w_v1 * w_v2 * (s_v - s) * exp(-w_v2 * pow(s_v - s, 2) - w_v4 * pow(d - d_v, 2));
				cout << d_v << "\t" << (s_v - s) << "\t" << 2 * w_v1 * w_v2 * (s_v - s) * exp(-w_v2 * pow(s - s_v, 2) - w_v4 * pow(d - d_v, 2)) << endl;
				//PF_grad_s += w_v1 * exp(- w_v4 * pow(d - d_v, 2)) / (s_v - s);
				//cout << (s - s_v) << "\t" << (d - d_v) << "\t" << w_v1 * exp(-w_v2 * pow(s - s_v, 2) - w_v4 * pow(d - d_v, 2)) * (s_v - s) / fabs(s_v - s) << endl;
				//PF_grad_d += w_v3 * exp(-w_v2 * pow(s - s_v, 2) - w_v4 * pow(d - d_v, 2)) * (d_v - d) / fabs(d_v - d);
				//if ((s - s_v) > -15 && (s - s_v) < 5)
				PF_grad_d += 2 * w_v1 * w_v4 * (d - d_v) * exp(-w_v2 * pow(s - s_v, 2) - w_v4 * pow(d - d_v, 2));
			}
		}

		//double PF_grad_s_norm = ref_velocity * PF_grad_s / sqrt(pow(PF_grad_s, 2) + pow(PF_grad_d, 2));
		//double PF_grad_d_norm = ref_velocity * PF_grad_d / sqrt(pow(PF_grad_s, 2) + pow(PF_grad_d, 2));
		if (PF_grad_s < -ref_velocity_s)
			PF_grad_s = -ref_velocity_s;
		if (PF_grad_s > 0)
			PF_grad_s = 0;
		if (fabs(PF_grad_d) > (5))
			PF_grad_d = 5 * PF_grad_d / fabs(PF_grad_d);
		
		PF_grad_d = 0;
		s -= PF_grad_s * 0.02;
		d -= PF_grad_d * 0.02;
		new_points = getXY(s, d, maps_s, maps_x, maps_y);

		xy_points.push_back({ s, d, -PF_grad_s, new_points[0], new_points[1] });

		speed_x = PF_grad_s;

		/*for (int j = 1; j < 11; j++) {
			s -= PF_grad_s * 0.02;
			d -= PF_grad_d * 0.02;
			cout << s << endl;
			new_points = getXY(s, d, maps_s, maps_x, maps_y);
			xy_points.push_back(new_points);
		}*/

		/*new_speed_x = (new_points[0] - x)/0.02;
		new_speed_y = (new_points[1] - y)/0.02;

		//limit speed
		speed_mag = sqrt(pow(new_speed_x, 2) + pow(new_speed_y, 2));
		if (speed_mag > ref_velocity_s) {
			new_speed_x = ref_velocity_s * new_speed_x / speed_mag;
			new_speed_y = ref_velocity_s * new_speed_y / speed_mag;
		}

		new_x = x + new_speed_x*0.02;
		new_y = y + new_speed_y*0.02;

		new_yaw = atan2(new_y - y, new_x - x);
		y = new_y;
		x = new_x;

		new_sd = getFrenet(x, y, new_yaw, maps_x, maps_y);
		s = new_sd[0]; 
		d = new_sd[1];*/
				
		//xy_points.push_back(new_points);
	}

	return xy_points;
}

double end_s, end_d, end_speed;

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

			//calculate ds and dd for each car in sensor_fusion
			vector<vector<double>> sensor_fusion_sd_frame;

			for (int i = 0; i < sensor_fusion.size() ; i++)
			{
				vector<double> sd_frame_components;

				double cur_s = sensor_fusion[i][5];
				double cur_d = sensor_fusion[i][6];

				double next_x = (double)sensor_fusion[i][1] + 0.02 * (double)sensor_fusion[i][3];
				double next_y = (double)sensor_fusion[i][2] + 0.02 * (double)sensor_fusion[i][4];

				auto next_sd = getFrenet(next_x, next_y, car_yaw, map_waypoints_x, map_waypoints_y);

				double ds = (next_sd[0] - cur_s) / 0.02;
				double dd = (next_sd[1] - cur_d) / 0.02;

				sd_frame_components.push_back(cur_s);
				sd_frame_components.push_back(cur_d);
				sd_frame_components.push_back(ds);
				sd_frame_components.push_back(dd);
				sensor_fusion_sd_frame.push_back(sd_frame_components);
			}

          	json msgJson;

			vector<double> next_x_vals;
			vector<double> next_y_vals;
			
			int previous_iteration = previous_path_x.size();
			//if (previous_iteration >= 0)
			//	previous_iteration = 0;

			if (previous_iteration < 2)
				previous_iteration = 0;

			//add previous path points
			for (int i = 0; i < previous_iteration; i++)
			{
				next_x_vals.push_back(previous_path_x[i]);
				next_y_vals.push_back(previous_path_y[i]);
			}

			if (previous_iteration >= 2)	{
				double prev_x = previous_path_x[previous_iteration - 1];
				double prev_xx = previous_path_x[previous_iteration - 2];
				double prev_y = previous_path_y[previous_iteration - 1];
				double prev_yy = previous_path_y[previous_iteration - 2];

				double prev_yaw = atan2(prev_y - prev_yy, prev_x - prev_xx);

				auto check_frenet = getFrenet(prev_x, prev_y, prev_yaw, map_waypoints_x, map_waypoints_y);
				end_path_s = check_frenet[0];
				end_path_d = check_frenet[1];
				//cout << "yaw" << "\t" << prev_yaw << endl;
				//cout << "sooo" << "\t" << end_path_s << "\t" << end_path_d << "\t" << prev_x << "\t" << prev_y << endl;
			} else {
				end_path_s = car_s;
				end_path_d = car_d;
				end_s = car_s;
				end_d = car_d;
				end_speed = car_speed;
			}

			auto new_path = potential_funtion(end_s, end_d, car_x, car_y, car_yaw, end_speed, 30, 100, previous_iteration, sensor_fusion_sd_frame, map_waypoints_s, map_waypoints_x, map_waypoints_y);

			end_s = new_path[new_path.size() - 1][0];
			end_d = new_path[new_path.size() - 1][1];
			end_speed = new_path[new_path.size() - 1][2];

			for (int i = 0; i < new_path.size(); i++)
			{
				next_x_vals.push_back(new_path[i][3]);
				next_y_vals.push_back(new_path[i][4]);
			}
			/*cout << "yeeeehaaaa" << endl;

			for (int i = 0; i < next_x_vals.size(); i++)
			{
				cout << next_x_vals[i] << "\t" << next_y_vals[i] << endl;
			}*/

			/*double next_s, next_d;

			double dist_increment = 0.5;
			for (int i = 0; i < 50; i++)
			{
				next_s = car_s + (i+1) * dist_increment;
				next_d = car_d;
				vector<double> next_xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
				next_x_vals.push_back(next_xy[0]);
				next_y_vals.push_back(next_xy[1]);
			}*/
          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;
          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}

