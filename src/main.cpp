#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

//maximum frenet s value of the highway
#define max_s 6945.554


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


//calculate distance between two points
double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}


//Find closest waypoint to specific point in map
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


//Find next waypoint to specific point and orientation in map
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


// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(vector<double> xvals, vector<double> yvals,	int order) {

	assert(xvals.size() == yvals.size());
	assert(order >= 1 && order <= xvals.size() - 1);

	Eigen::MatrixXd A(xvals.size(), order + 1);

	for (unsigned int i = 0; i < xvals.size(); i++) {
		A(i, 0) = 1.0;
	}

	for (unsigned int j = 0; j < xvals.size(); j++) {
		for (unsigned int i = 0; i < order; i++) {
			A(j, i + 1) = A(j, i) * xvals[j];
		}
	}

	auto Q = A.householderQr();
	Eigen::VectorXd Yvals = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(yvals.data(), yvals.size());
	auto result = Q.solve(Yvals);
	return result;
}


// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
	double result = 0.0;
	for (unsigned int i = 0; i < coeffs.size(); i++) {
		result += coeffs[i] * pow(x, i);
	}
	return result;
}


// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double dx, double dy, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

	int prev_wp;
	prev_wp = next_wp - 1;
	if (next_wp == 0)
	{
		prev_wp = maps_x.size() - 1;
	}

	double n_x = maps_x[next_wp] - maps_x[prev_wp];
	double n_y = maps_y[next_wp] - maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x + x_y * n_y) / (n_x*n_x + n_y * n_y);
	double proj_x = proj_norm * n_x;
	double proj_y = proj_norm * n_y;

	//distance from projection
	double frenet_d = distance(x_x, x_y, proj_x, proj_y);

	//find the projection of dx onto n
	double dproj_norm = (dx*n_x + dy * n_y) / (n_x*n_x + n_y * n_y);
	double dproj_x = dproj_norm * n_x;
	double dproj_y = dproj_norm * n_y;

	//rate of change of distance from projection
	double frenet_dd = distance(dx, dy, dproj_x, dproj_y);

	//see if d value is positive or negative by comparing it to a center point
	double center_x = 1000 - maps_x[prev_wp];
	double center_y = 2000 - maps_y[prev_wp];
	double centerToPos = distance(center_x, center_y, x_x, x_y);
	double centerToRef = distance(center_x, center_y, proj_x, proj_y);

	if (centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	//see if dd value is positive or negative by comparing it to a center point
	double dcenterToPos = distance(center_x, center_y, dx, dy);
	double dcenterToRef = distance(center_x, center_y, dproj_x, dproj_y);

	if (dcenterToPos <= dcenterToRef)
	{
		frenet_dd *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for (int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
	}

	frenet_s += distance(0, 0, proj_x, proj_y);

	//calculate ds value
	double frenet_ds = distance(0, 0, dproj_x, dproj_y);

	return { frenet_s, frenet_d, frenet_ds, frenet_dd };
}


// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
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

	return {x,y};
}


//calculate cost of driving at each of the three lanes of the highway
vector<double> lane_cost(double start_s, double start_d, double ref_velocity, double cur_velocity, vector<vector<double>> sensor_fusion_sd_frame)
{
	vector<double> costs;

	//Cost function of driving at each lane
	//CF =   time integration from 0 to 15 ( summation of all vehciles ( w1 * exp( - w2 * delta_s^2 - w3 * delta_d^2) ) exp(- w4 * t^2) )
	//This cost function penalizes driving near vehicles as time progresses

	//initializing weights in the cost function
	double w1 = 10;
	double w2 = 0.05;
	double w3 = 0.3;
	double w4 = 0.003;
	

	//interate through lanes
	for (int i = 0 ; i < 3 ; i++)
	{

		double d = 4 * i + 2; //frenet d value of driving at the center of the lane

		double cost = 0; //initialize cost
		double s = start_s; //initial frenet s value
		double car_velocity = 0.44704 * cur_velocity; //initial car velocity (converted from mph to m/s)
		
		//integrate over time as vehicles movement progresses
		for (int t = 0; t < 15 ; t++)
		{
			//reset s value if it exceeds max value
			while (s > max_s)
				s -= max_s;

			//iterate through all detected vehicles
			for (int j = 0; j < sensor_fusion_sd_frame.size(); j++)
			{
				
				//initial values of other vehicles in frenet coordinates
				double start_s_v = sensor_fusion_sd_frame[j][0]; 
				double start_d_v = sensor_fusion_sd_frame[j][1];
				
				//ignoring car on the other side of the road
				if (start_d_v > 0)
				{
					//ignoring cars behind our vehicle in the same lane
					if ((start_s_v > start_s) || (fabs(start_d_v - start_d) > 1)) {

						//linearly interpolate the position of each car
						double s_v = sensor_fusion_sd_frame[j][0] + t * sensor_fusion_sd_frame[j][2];
						double d_v = sensor_fusion_sd_frame[j][1] + t * sensor_fusion_sd_frame[j][3];

						while (s_v > max_s)
							s_v -= max_s;

						//update cost function
						cost += w1 * exp(-w2 * pow(s - s_v, 2) - w3 * pow(d - d_v, 2)) * exp(- w4 * t * t);
					}
				}
			
			}
			
			//update position of car assuming its accelerating towards reference velocity
			s += car_velocity;
			if (car_velocity < ref_velocity)
			{
				car_velocity += 5;
				s += 0.5 * 5;
			}
		}

		costs.push_back(cost);
	}
	
	//return vector of cost of each lane
	return costs;
}


//determine best lane and frenet d value to drive on based on cost of each lane
double desired_d(double car_d, vector<double> lane_cost)
{
	//properties of current lane
	int car_lane = (int)car_d / 4; 
	double car_lane_cost = lane_cost[car_lane];

	//find minimum cost lane
	int desired_lane;
	double min_lane_cost = std::numeric_limits<double>::infinity();

	for (int i = 0; i < lane_cost.size(); i++)
	{
		if (lane_cost[i] < min_lane_cost)
		{
			min_lane_cost = lane_cost[i];
			desired_lane = i;
		}
	}

	//in case double lane change is needed, transfer to mid lane first
	if (fabs(desired_lane - car_lane) > 1)
	{
		desired_lane = car_lane + (desired_lane - car_lane) / fabs(desired_lane - car_lane);
		min_lane_cost = lane_cost[desired_lane];
	}

	//leaving a deadband (as value and precentage) for current lane and ensure lane transfer is safe
	if ((min_lane_cost > (0.5 * car_lane_cost)) || ( (car_lane_cost - min_lane_cost) < 4 ) || (min_lane_cost > 10))
		desired_lane = car_lane;	

	//return desired frenet d value
	return 4 * desired_lane + 2;
}



//generate waypoints based on current frenet position and desired lane and velocity
vector<vector<double>> waypoints(double s, double d, double desired_d, double ref_velocity, double car_velocity, vector<vector<double>> sensor_fusion_sd_frame, int start_point, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	vector<vector<double>> xy_points; //vector for x and y points of the desired trajectory
	
	//start from previous path end point to reach 50 way points
	for (int i = start_point; i < 50; i++)
	{
		double acc_status = 0; //acceleration status: -1--> deaccelerate, 0--> keep speed, 1--> accelerate

		if (car_velocity < ref_velocity)
			acc_status = 1;

		//iterate through cars to ensure safe waypoints
		for (int j = 0; j < sensor_fusion_sd_frame.size(); j++)
		{

			//linearly interpolate the position of each car
			double s_v = sensor_fusion_sd_frame[j][0] + i * 0.02 * sensor_fusion_sd_frame[j][2];
			double d_v = sensor_fusion_sd_frame[j][1] + i * 0.02 * sensor_fusion_sd_frame[j][3];

			//limit frenet s value
			while (s_v > max_s)
				s_v -= max_s;

			//determine acceleration status based on vehicles around
			if ((s_v > s) && ((s_v - s) < 25) && (fabs(d_v - d) < 2))
			{
				acc_status = 0;
 				if ( ((s_v - s) < 20) && (car_velocity > sensor_fusion_sd_frame[j][2]))
					acc_status = -1;
			}
		}	

		//accelerate car in s direction
		car_velocity += (acc_status * 0.1);
		if (car_velocity < 0)
			car_velocity = 0;
		
		//up waypoints s and d values
		s += car_velocity * 0.02;
		d += 0.02 * (desired_d - d) * car_velocity / 20;
		
		while (s > max_s)
			s -= max_s;

		//convert to cartessian coordinates and add to waypoint vector
		auto xy_point = getXY(s, d, maps_s, maps_x, maps_y);

		xy_points.push_back({ s, d, car_velocity, xy_point[0], xy_point[1] });
	}

	return xy_points;
}

//convert sensor fusion data from cartessian to frenet to find ds and dd
vector<vector<double>> convert_sensor_fusion(vector<vector<double>> sensor_fusion, const vector<double> &maps_x, const vector<double> &maps_y)
{
	vector<vector<double>> sensor_fusion_sd_frame;

	//iterate through detected vehicles
	for (int i = 0; i < sensor_fusion.size(); i++)
	{
		vector<double> sd_frame_components;

		double cur_x = sensor_fusion[i][1];
		double cur_y = sensor_fusion[i][2];
		double cur_dx = sensor_fusion[i][3];
		double cur_dy = sensor_fusion[i][4];
		double cur_s = sensor_fusion[i][5];
		double cur_d = sensor_fusion[i][6];

		double cur_yaw = atan2(cur_dy, cur_dx);

		auto cur_sd = getFrenet(cur_x, cur_y, cur_dx, cur_dy, cur_yaw, maps_x, maps_y);

		double ds = cur_sd[2];
		double dd = cur_sd[3];
		
		sd_frame_components.push_back(cur_s);
		sd_frame_components.push_back(cur_d);
		sd_frame_components.push_back(ds);
		sd_frame_components.push_back(dd);

		sensor_fusion_sd_frame.push_back(sd_frame_components);
	}

	return sensor_fusion_sd_frame;
}

//global variables
double end_s, end_d, end_speed; //final frenet coordinates and speed of last point on previous path
double d_target; //target frenet d value
int iteration = 0; //number of iteration
vector<double> accum_lane_cost = { 0, 0, 0 }; //accumulated lane cost
double ref_velocity = 22; //reference velocity

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


			//convert sensor fusion data to frenet coordinates (to find ds and dt)
			auto sensor_fusion_sd_frame = convert_sensor_fusion(sensor_fusion, map_waypoints_x, map_waypoints_y);

			json msgJson;

			vector<double> next_x, next_y; //x and y path waypoints 
			vector<double> next_x_vals, next_y_vals; //final trajectory points after processing

			//use a maximum of 25 points from previous path
			double last_point = previous_path_x.size();
			if (last_point > 25)
				last_point = 25;

			//add remaining path points to current trajectory
			for (int i = 0; i < previous_path_x.size(); i++)
			{
				//points for polynomial fitting
				next_x.push_back(previous_path_x[i]);
				next_y.push_back(previous_path_y[i]);
				if (i < last_point) {
					//directly feed some points to next trajectoy
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
				}
			}

			//initialize previous path end points in frenet coordinates
			if (previous_path_x.size() == 0)
			{
				end_s = car_s;
				end_d = car_d;
				end_speed = 0.44704*car_speed;
			}

			//calculate cost of each lane
			auto cost = lane_cost(car_s, car_d, ref_velocity, car_speed, sensor_fusion_sd_frame);

			//accumulate newly calculated cost
			for (int i = 0; i < cost.size(); i++)
				accum_lane_cost[i] += cost[i];

			if ((iteration % 100 == 0))
			{
				//average the accumulated costs
				for (int i = 0; i < cost.size(); i++)
					accum_lane_cost[i] = accum_lane_cost[i] / 100.0;

				//find best frenet d value for driving
				d_target = desired_d(end_d, accum_lane_cost);
				
				//reset accumulated cost values
				accum_lane_cost = { 0, 0, 0 };
			}
			
			//update waypoints
			auto path = waypoints(end_s, end_d, d_target, ref_velocity, end_speed, sensor_fusion_sd_frame, previous_path_x.size(), map_waypoints_s, map_waypoints_x, map_waypoints_y);
			
			//update previous path end points
			if (path.size()) {
				end_s = path[path.size() - 1][0];
				end_d = path[path.size() - 1][1];
				end_speed = path[path.size() - 1][2];
			}
			
			//concatenate newly calculated 
			for (int i = 0; i < path.size(); i++)
			{
				next_x.push_back(path[i][3]);
				next_y.push_back(path[i][4]);
			}

			//select few waypoints for time polynomial fitting
			vector<double> ptsx, ptsy; //selected waypoints for polynomial fitting
			vector<double> time;
			for (int i = 0; i < next_x.size(); i++)
			{
				if ( ((i % 25) == 0) || (i==(next_x.size()-1)) )
				{
					ptsx.push_back(next_x[i]);
					ptsy.push_back(next_y[i]);
					time.push_back((double)(i+1.0) * 0.02);
				}
			}

			//fit selected points to a second order polynomial
			auto x_poly = polyfit(time, ptsx, 2);
			auto y_poly = polyfit(time, ptsy, 2);
			
			//add points from fitted polynomial to final trajectory points
			for (int i = (1+last_point); i < 51; i++)
			{
				//evalute x and y points
				double t = (double)i * 0.02;
				double xcar = polyeval(x_poly, t);
				double ycar = polyeval(y_poly, t);
				
				//limit velocity not to exceed limit

				double vel_x, vel_y; //velocity in x and y directions

				//calculate velocity
				if (next_x_vals.size() > 1) {
					vel_x = (xcar - next_x_vals[next_x_vals.size() - 1])/0.02;
					vel_y = (ycar - next_y_vals[next_y_vals.size() - 1])/0.02;
				}
				else {
					vel_x = (xcar - car_x)/0.02;
					vel_y = (ycar - car_y)/0.02;
				}
				
				//limit velocity
				double speed = sqrt(pow(vel_x, 2) + pow(vel_y, 2));
				if (speed > ref_velocity) {
					vel_x = vel_x * ref_velocity / speed;
					vel_y = vel_y * ref_velocity / speed;
				}
				
				//update points based on limited velocity
				if (next_x_vals.size() > 1) {
					xcar = next_x_vals[next_x_vals.size() - 1] + vel_x * 0.02;
					ycar = next_y_vals[next_y_vals.size() - 1] + vel_y * 0.02;
				}
				else {
					xcar = car_x + vel_x * 0.02;
					ycar = car_y + vel_y * 0.02;
				}

				//push points
				next_x_vals.push_back(xcar);
				next_y_vals.push_back(ycar);

			}

			//reevaluate end path points to mitigate large biases 
			//due to continuous conversion between cartessian and frenet
			if ((iteration % 500) == 0)
			{
				double xdiff = next_x_vals[next_x_vals.size() - 1] - next_x_vals[next_x_vals.size() - 2];
				double ydiff = next_y_vals[next_y_vals.size() - 1] - next_y_vals[next_y_vals.size() - 2];
				double final_yaw = atan2(ydiff, xdiff);
				auto end_frenet = getFrenet(next_x_vals[next_x_vals.size() - 1], next_y_vals[next_y_vals.size() - 1], xdiff / 0.02, ydiff / 0.02, final_yaw, map_waypoints_x, map_waypoints_y);
				end_s = end_frenet[0];
				end_d = end_frenet[1];
				end_speed = end_frenet[2];
			}

			iteration++; //increment number of iterations
			
			//pass waypoint to Json message
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

