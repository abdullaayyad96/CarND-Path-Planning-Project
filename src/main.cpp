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

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenetss(double x, double y, double dx, double dy, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
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

	double dproj_norm = (dx*n_x + dy * n_y) / (n_x*n_x + n_y * n_y);
	double dproj_x = dproj_norm * n_x;
	double dproj_y = dproj_norm * n_y;

	double frenet_d = distance(x_x, x_y, proj_x, proj_y);
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


vector<double> lane_cost(double start_s, double start_d, double ref_velocity, double cur_velocity, vector<vector<double>> sensor_fusion_sd_frame)
{
	//return cost of driving at each lane
	vector<double> costs;

	//Construct cost function
	//CF =   summation of w_v1 exp(- w_v2 * (s - s_v)^2  - w_v3 * (d - d_v)^2 ) for every vehicle

	double w_v1 = 10;
	double w_v2 = 0.05;
	double w_v3 = 0.3;
	

	for (int i = 0 ; i < 3 ; i++)
	{
		//iterate through lanes
		double d = 4 * i + 2;
		double cost = 0;
		double s = start_s;
		double car_velocity = 0.44704 * cur_velocity;

		for (int t = 0; t < 15 ; t++)
		{
			//integrate for 10 seconds 

			s = start_s + t * car_velocity;

			while (s > 6945.554)
				s -= 6945.554;

			if (car_velocity < ref_velocity)
				car_velocity += 0.1;
			
			for (int j = 0; j < sensor_fusion_sd_frame.size(); j++)
			{
				//iterate through cars
				double start_s_v = sensor_fusion_sd_frame[j][0];
				double d_v = sensor_fusion_sd_frame[j][1];// +t * sensor_fusion_sd_frame[j][3];

				if ((start_s_v > start_s) || (fabs(d_v - start_d) > 1)) {

					//linearly interpolate the position of each car
					double s_v = sensor_fusion_sd_frame[j][0] + t * sensor_fusion_sd_frame[j][2];

					while (s_v > 6945.554)
						s_v -= 6945.554;

					//if (fabs(s_v - s) < 5)
					//	cout << j << "\t" << t << "\t" << (s_v - s) << "\t" << d_v << endl;

					//update PF
					if (d_v > 0)
						cost += w_v1 * exp(-w_v2 * pow(s - s_v, 2) - w_v3 * pow(d - d_v, 2)) *exp(-0.005 * t * t);
				}
			
			}
		}
		//cout << d << "\t" << cost << endl;
		costs.push_back(cost);
		
	}
	return costs;
}


double desired_d(double car_d, vector<double> lane_cost)
{
	//return desired d value for a lane

	int car_lane = (int)car_d / 4;
	double car_lane_cost = lane_cost[car_lane];

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

	//leaving a deadband for current lane
	if ((min_lane_cost > (0.7 * car_lane_cost)) || ( (car_lane_cost - min_lane_cost) < 5 ))
		desired_lane = car_lane;	

	//in case double lane change is needed, make sure mid lane is safe
	if ((fabs(desired_lane - car_lane) > 1) && (lane_cost[1] > 0.4 * lane_cost[car_lane]))
		desired_lane = car_lane;
	else if ((fabs(desired_lane - car_lane) > 1) && ( (min_lane_cost > 0.7 * lane_cost[1]) || ((lane_cost[1] - min_lane_cost) < 5)) )
		desired_lane = 1;


	return 4 * desired_lane + 2;
}

vector<vector<double>> trajectory(double s, double d, double desired_d, double ref_velocity, double car_velocity, vector<vector<double>> sensor_fusion_sd_frame, int start_point, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	vector<vector<double>> xy_points;

	double w_v1 = 100;
	double w_v2 = 0.05;
	double w_v3 = 1;

	for (int i = start_point; i < 50; i++)
	{
		if (car_velocity < ref_velocity)
		{
			car_velocity += 0.1;
		}

		for (int j = 0; j < sensor_fusion_sd_frame.size(); j++)
		{
			//iterate through cars

			//linearly interpolate the position of each car
			double s_v = sensor_fusion_sd_frame[j][0] + i * 0.02 * sensor_fusion_sd_frame[j][2];
			double d_v = sensor_fusion_sd_frame[j][1];// +i * 0.02 * sensor_fusion_sd_frame[j][3];

			while (s_v > 6945.554)
				s_v -= 6945.554;

			//update PF
			//if ( (s_v > s) && ( fabs(d_v-d) < 2 ) )
			//	car_velocity += -1.5 * exp(-0.005 * pow(s_v - s, 2));// *exp(-10 * pow(d_v - d, 2));
			if ( (s_v > s) && ((s_v - s) < 15) && (fabs(d_v - d) < 2))
				car_velocity += -1.5 + 1.5*(s_v - s)/15;
		}	

		if (car_velocity < 0)
			car_velocity = 0;

		/*if ( fabs(desired_d - d) > 1 )
			d += 1 * 0.02 * (desired_d - d) / fabs(desired_d - d);
		else 
			d += 0.02 * (desired_d - d);*/

		d += 0.02 * (desired_d - d) * car_velocity / 20;

		s += car_velocity * 0.02;

		while (s > 6945.554)
			s -= 6945.554;

		auto xy_point = getXY(s, d, maps_s, maps_x, maps_y);

		xy_points.push_back({ s, d, car_velocity, xy_point[0], xy_point[1] });
	}

	return xy_points;
}

double end_s, end_d, end_speed, d_target;
int iteration = 0;
vector<double> speeds;
vector<double> accum_lane_cost = { 0, 0, 0 };
double startacc_x = 0;
double startacc_y = 0;
double ref_velocity = 20;
double max_acc = 10;
vector<double> x_poly, y_poly;
double last_vel_x, last_vel_y, prev_vel_x = 0, prev_vel_y = 0, last_acc_x, last_acc_y, last_speed, last_acc;


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

			for (int i = 0; i < sensor_fusion.size(); i++)
			{
				vector<double> sd_frame_components;

				double cur_s = sensor_fusion[i][5];
				double cur_d = sensor_fusion[i][6];

				double next_x = (double)sensor_fusion[i][1] + 0.02 * (double)sensor_fusion[i][3];
				double next_y = (double)sensor_fusion[i][2] + 0.02 * (double)sensor_fusion[i][4];

				auto next_sd = getFrenet(next_x, next_y, deg2rad(car_yaw), map_waypoints_x, map_waypoints_y);

				double ds = (next_sd[0] - cur_s) / 0.02;
				double dd = (next_sd[1] - cur_d) / 0.02;

				/*cout << "ok" << endl;
				cout << cur_d << "\t" << next_sd[1] << endl;
				cout << "ok" << endl;*/
				double cur_speed = sqrt( pow((double)sensor_fusion[i][3], 2) + pow((double)sensor_fusion[i][4], 2) );

				sd_frame_components.push_back(cur_s);
				sd_frame_components.push_back(cur_d);
				sd_frame_components.push_back(cur_speed);
				sd_frame_components.push_back(ds);
				sd_frame_components.push_back(dd);
				sensor_fusion_sd_frame.push_back(sd_frame_components);
			}

			json msgJson;

			vector<double> ptsx;
			vector<double> ptsy;
			vector<double> next_x;
			vector<double> next_y;
			vector<double> next_x_vals;
			vector<double> next_y_vals;
			
			double last_point = previous_path_x.size();
			if (last_point > 25)
				last_point = 25;

			for (int i = 0; i < previous_path_x.size(); i++)
			{
				next_x.push_back(previous_path_x[i]);
				next_y.push_back(previous_path_y[i]);
				if (i < 25) {
					next_x_vals.push_back(previous_path_x[i]);
					next_y_vals.push_back(previous_path_y[i]);
				}
			}


			if (previous_path_x.size() == 0) 
			{
				end_path_s = car_s;
				end_path_d = car_d;
				end_s = car_s;
				end_d = car_d;
				end_speed = 0.44704*car_speed;
			}

			auto cost = lane_cost(car_s, car_d, 20, car_speed, sensor_fusion_sd_frame);
			for (int i = 0; i < cost.size(); i++)
				accum_lane_cost[i] += cost[i];
			if ((iteration % 100 == 0))
			{
				for (int i = 0; i < cost.size(); i++)
					accum_lane_cost[i] = accum_lane_cost[i] / 100.0;

				d_target = desired_d(end_d, accum_lane_cost);

				cout << accum_lane_cost[0] << "\t" << accum_lane_cost[1] << "\t" << accum_lane_cost[2] << endl;
				cout << d_target << endl;
				accum_lane_cost = { 0, 0, 0 };
			}
			
			auto path = trajectory(end_s, end_d, d_target, ref_velocity, end_speed, sensor_fusion_sd_frame, previous_path_x.size(), map_waypoints_s, map_waypoints_x, map_waypoints_y);

			iteration++;
			/*if (path.size()) {
				end_s = path[path.size() - 1][0];
				end_d = path[path.size() - 1][1];
				end_speed = path[path.size() - 1][2];
				//end_speed = sqrt(prev_vel_x * prev_vel_x + prev_vel_y * prev_vel_y);
				//cout << "speeds \t" << prev_vel_x << "\t" << prev_vel_y << "\t" << end_speed << endl;
			}*/

			for (int i = 0; i < path.size(); i++)
			{
				next_x.push_back(path[i][3]);
				next_y.push_back(path[i][4]);
				speeds.push_back(path[i][2]);
				if (speeds.size() > 50)
					speeds.erase(speeds.begin());
			}

			vector<double> time;
			for (int i = 0; i < next_x.size(); i++)
			{
				if ( ((i % 25) == 0) || (i==(next_x.size()-1)) )
				{
					ptsx.push_back(next_x[i]);
					ptsy.push_back(next_y[i]);
					//cout << "fit \t" << next_x[i] << "\t" << next_y[i] << endl;
					//ptsx.push_back((next_x[i] - car_x)*cos(deg2rad(car_yaw)) + (next_y[i] - car_y)*sin(deg2rad(car_yaw)));
					//ptsy.push_back(-(next_x[i] - car_x)*sin(deg2rad(car_yaw)) + (next_y[i] - car_y)*cos(deg2rad(car_yaw)));
					time.push_back((double)(i+1.0) * 0.02);
				}
			}

			auto x_poly = polyfit(time, ptsx, 2);
			auto y_poly = polyfit(time, ptsy, 2);
			/*auto xypoly = polyfit(ptsx, ptsy, 3);

			double start_x;
			if (previous_path_x.size() == 0)
				start_x = 0;
			else
				start_x = ((double)previous_path_x[previous_path_x.size()-1] - car_x)*cos(deg2rad(car_yaw)) + ((double)previous_path_y[previous_path_y.size() - 1] - car_y)*sin(deg2rad(car_yaw));

			for (int i = (last_point); i < 5 ; i++)
			{
				double xcar = start_x + 0.02 * speeds[i];
				double ycar = polyeval(xypoly, xcar);
				//cout << xcar << "\t" << ycar << endl;
				start_x = xcar;

				next_x_vals.push_back(car_x + xcar * cos(deg2rad(car_yaw)) - ycar * sin(deg2rad(car_yaw)));
				next_y_vals.push_back(car_y + ycar * cos(deg2rad(car_yaw)) + xcar * sin(deg2rad(car_yaw)));
			}
			*/


			//third order time polynomials / minimum jerk polynomials 
			
			//if (previous_path_x.size() != 0)
			//{
			//	startacc_x = JMT_eval_acc(x_poly, 1 - 0.02*previous_path_x.size());
			//	startacc_y = JMT_eval_acc(y_poly, 1 - 0.02*previous_path_x.size());
			//}
			//double final_yaw = atan2(next_y[next_x.size() - 1] - next_y[next_x.size() - 2], next_x[next_x.size() - 1] - next_x[next_x.size() - 2]);

			//vector<double> xstart = { car_x, 0.44704*car_speed*cos(deg2rad(car_yaw)), startacc_x };
			//vector<double> xend = { next_x[next_x.size() - 1], end_speed*cos(final_yaw), 0};
			//x_poly = JMT(xstart, xend, 0.02*(next_x.size()));
			//vector<double> ystart = { car_y, 0.44704*car_speed*sin(deg2rad(car_yaw)), startacc_y };
			//vector<double> yend = { next_y[next_y.size() - 1], end_speed*sin(final_yaw), 0 };
			//y_poly = JMT(ystart, yend, 0.02*(next_x.size()));
			double xcar, ycar;
			//cout << "speeds" << car_speed * cos(deg2rad(car_yaw)) << "\t" << end_speed * cos(final_yaw) << "\t" << car_speed * cos(deg2rad(car_yaw)) << "\t" << end_speed * sin(final_yaw) << endl;

			for (int i = (1+last_point); i < 51; i++)
			{
				double t = (double)i * 0.02;
				xcar = polyeval(x_poly, t);
				ycar = polyeval(y_poly, t);
				
				if (next_x_vals.size() > 1) {
					last_vel_x = (xcar - next_x_vals[next_x_vals.size() - 1])/0.02;
					last_vel_y = (ycar - next_y_vals[next_y_vals.size() - 1])/0.02;
				}
				else {
					last_vel_x = (xcar - car_x)/0.02;
					last_vel_y = (ycar - car_y)/0.02;
				}
				
				//cout << "1 \t" << last_vel_x << "\t" << last_vel_y << endl;
				last_speed = sqrt(pow(last_vel_x, 2) + pow(last_vel_y, 2));
				//cout << last_speed << endl;
				if (last_speed > ref_velocity) {
					last_vel_x = last_vel_x * ref_velocity / last_speed;
					last_vel_y = last_vel_y * ref_velocity / last_speed;
				}
				//cout << "2 \t" << last_vel_x << "\t" << last_vel_y << endl;
				/*
				last_acc_x = (last_vel_x - prev_vel_x) / 0.02;
				last_acc_y = (last_vel_y - prev_vel_y) / 0.02;
				last_acc = sqrt(pow(last_acc_x, 2) + pow(last_acc_y, 2));
				//cout << "1 \t" << last_vel_x << "\t" << last_vel_y << endl;
				if (last_acc > max_acc) {
					last_acc_x = last_acc_x * max_acc / last_acc;
					last_acc_y = last_acc_y * max_acc / last_acc;
					last_vel_x = prev_vel_x + 0.02 * last_acc_x;
					last_vel_y = prev_vel_y + 0.02 * last_acc_y; 
				}
				/
				//cout << "2 \t" << last_vel_x << "\t" << last_vel_y << endl;
				*/
				if (next_x_vals.size() > 1) {
					xcar = next_x_vals[next_x_vals.size() - 1] + last_vel_x * 0.02;
					ycar = next_y_vals[next_y_vals.size() - 1] + last_vel_y * 0.02;
				}
				else {
					xcar = car_x + last_vel_x * 0.02;
					ycar = car_y + last_vel_y * 0.02;
				}
				prev_vel_x = last_vel_x;
				prev_vel_y = last_vel_y;
				
					
				//xcar = JMT_eval(x_poly, t);
				//ycar = JMT_eval(y_poly, t);

				next_x_vals.push_back(xcar);
				next_y_vals.push_back(ycar);

			}
			
			double xdiff = next_x_vals[next_x_vals.size() - 1] - next_x_vals[next_x_vals.size() - 2];
			double ydiff = next_y_vals[next_y_vals.size() - 1] - next_y_vals[next_y_vals.size() - 2];
			double final_yaw = atan2(ydiff, xdiff);
			auto end_frenet = getFrenetss(next_x_vals[next_x_vals.size() - 1], next_y_vals[next_y_vals.size() - 1], xdiff / 0.02, ydiff / 0.02, final_yaw, map_waypoints_x, map_waypoints_y);
			end_s = end_frenet[0];
			end_d = end_frenet[1];
			if (next_x_vals.size() > previous_path_x.size())
				end_speed = path[next_x_vals.size() - previous_path_x.size() - 1][2];

			
			
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

