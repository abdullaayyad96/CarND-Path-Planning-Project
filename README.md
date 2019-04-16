# CarND-Path-Planning-Project
Self-Driving Car Engineer Nanodegree Program


[//]: # (Picture Definition)
[cost]: https://user-images.githubusercontent.com/37302013/48884687-e02f3600-ee68-11e8-89f4-5d2d9961a411.png
[frenet]: https://user-images.githubusercontent.com/37302013/48889248-8c2c4d80-ee78-11e8-9d3d-71ba9f0afaea.png
[test]: https://user-images.githubusercontent.com/37302013/48889855-80418b00-ee7a-11e8-8401-3e3b57b61acb.png

## Project Discription

This project is part of UDacity's Self-Driving Car Engineer Nanodegree Program. The purpose is to develop a path/trajectory planner to safely navigate a vehicle in a highway.

### Simulator.
You can download the Term3 Simulator which contains the Path Planning Project from the [releases tab (https://github.com/udacity/self-driving-car-sim/releases/tag/T3_v1.2).

### Goals
In this project, the goal is to safely navigate around a virtual highway with other traffic that is driving +-10 MPH of the 50 MPH speed limit. The path planner is provided with the car's localization and sensor fusion data, there is also a sparse map list of waypoints around the highway. The car should try to go as close as possible to the 50 MPH speed limit, which means passing slower traffic when possible, note that other cars will try to change lanes too. The car should avoid hitting other cars at all cost as well as driving inside of the marked road lanes at all times, unless going from one lane to another. The car should be able to make one complete loop around the 6946m highway. Since the car is trying to go 50 MPH, it should take a little over 5 minutes to complete 1 loop. Also the car should not experience total acceleration over 10 m/s^2 and jerk that is greater than 10 m/s^3.

#### The map of the highway is in data/highway_map.txt
Each waypoint in the list contains  [x,y,s,dx,dy] values. x and y are the waypoint's map coordinate position, the s value is the distance along the road to get to that waypoint in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

The highway's waypoints loop around so the frenet s value, distance along the road, goes from 0 to 6945.554.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

Here is the data provided from the Simulator to the C++ Path Planner Program

#### Main car's localization Data (No Noise)

["x"] The car's x position in map coordinates

["y"] The car's y position in map coordinates

["s"] The car's s position in frenet coordinates

["d"] The car's d position in frenet coordinates

["yaw"] The car's yaw angle in the map

["speed"] The car's speed in MPH

#### Previous path data given to the Planner

//Note: Return the previous list but with processed points removed, can be a nice tool to show how far along
the path has processed since last time. 

["previous_path_x"] The previous list of x points previously given to the simulator

["previous_path_y"] The previous list of y points previously given to the simulator

#### Previous path's end s and d values 

["end_path_s"] The previous list's last point's frenet s value

["end_path_d"] The previous list's last point's frenet d value

#### Sensor Fusion Data, a list of all other car's attributes on the same side of the road. (No Noise)

["sensor_fusion"] A 2d vector of cars and then that car's [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates. 

## Details

1. The car uses a perfect controller and will visit every (x,y) point it recieves in the list every .02 seconds. The units for the (x,y) points are in meters and the spacing of the points determines the speed of the car. The vector going from a point to the next point in the list dictates the angle of the car. Acceleration both in the tangential and normal directions is measured along with the jerk, the rate of change of total Acceleration. The (x,y) point paths that the planner recieves should not have a total acceleration that goes over 10 m/s^2, also the jerk should not go over 50 m/s^3. (NOTE: As this is BETA, these requirements might change. Also currently jerk is over a .02 second interval, it would probably be better to average total acceleration over 1 second and measure jerk from that.

2. There will be some latency between the simulator running and the path planner returning a path, with optimized code usually its not very long maybe just 1-3 time steps. During this delay the simulator will continue using points that it was last given, because of this its a good idea to store the last points you have used so you can have a smooth transition. previous_path_x, and previous_path_y can be helpful for this transition since they show the last points given to the simulator controller with the processed points already removed. You would either return a path that extends this previous path or make sure to create a new path that has a smooth transition with this last path.

## Dependencies

* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

## Project Instructions and Rubric

The project has been developed to follow this [rubric](https://review.udacity.com/#!/rubrics/1971/view).

## Trajectory Generation Methodology

While driving in a highway might seem simple at first glance, a deeper look unfolds the complexity of the problem. The main challenges manifest due to the existence of the other vehicles that can obstruct the vehicle's progress raising the need for a prediction and behavior planning steps to determine the most appropriate lane and speed to drive at. Once these are determined, a safe and executable trajectory should be determined for the vehicle to follow.

My approach divides these process into four main steps discussed below:

### Lane Cost Calculation:

Determining the best lane to drive at first requires a criteria for numerically evaluating the appropriateness of each lane. For this purpose, an exponential cost function has been developed as:

![cost_function][cost]

s: target vehicle frenet s value
d: target vehicle frenet d value
s_v: other vehicles in the road frenet s value
d_v: other vehicles in the road frenet d value
w1, w2, w3, w4: weights

(for brief description of frenet coordinates, check final section)

The cost function clearly penalizes driving near other vehicles in the road by incorporating the difference in the frenet coordinates. In order to calculate the cost of driving at each lane, the calculation is repeated several times with varying the d value to represent driving at the center of each lane.

A time horizon was also implemented as a way of allowing a form of future prediction. This means that the position of all vehicles need to be updated as time progresses. This is done by linearly interpolating the position of other vehicles based on their last measured velocity. While for the target vehicle, the position is updated as if the vehicle is accelerating towards the target speed. While this prediction is straightforward and does not predict complex behaviors of other vehicles, it provides very good results as it's iteratively applied with a short time constant.

Finally, another exponential term related to time is added so that higher priority is given to the closest time steps as opposed to the furthest time steps which can be compensated for by future predictions and behavior planning.

This sub-process is implemented in the lane_cost function in lines 235-310 of "srs/main.cpp".

### Determining desired lane

The second step of the path planning algorithm is deciding which path is best to drive at which is the most straightforward step of the algorithm. First, the lane costs calculated from the previous steps are averaged over the number of times they have been calculated. The averaged lane costs are then used to determine which lane is most appropriate. This included ensuring that lane changes are safe, determining the feasibility and worthiness of lane changes and ensuring that lane changes are applied one lane at a time. This is implemented in desired_d function in lines 313-345 of "srs/main.cpp".

### Waypoint generation

Once a desired lane is determined, a rough trajectory is generated in order to move the car to that lane while accelerating to the desired speed. The generated waypoints also need to ensure that car does not collide with other vehicle on the road. Initially, these waypoints are generated in frenet coordinates and then converted to cartessian coordinates. It's worth noting that at every iteration, some remaining waypoints from the last generated path are taken into account and the new waypoints are simply aggregated. Most of this part is implemented in the waypoints function in lines 350-402 of "srs/main.cpp". While the waypoints transfer from previous path are impelemented in lines 531-543.

### Smoothening the waypoints with polynomial fitting

The last step is smoothening the waypoints. The main reason why smoothening is needed is because in the previous steps, waypoint generation is first done in frenet and then converted to cartesian coordinates. This leads to the generated path being indifferentiable at some points. This causes very high and unrealistic accelerations at times. Therefore, several waypoints are selected and fitted into second order polynomials ensuring a continuously differentiable trajectory/path. Waypoints are then sampled from these polynomials with a 0.02s time difference and fed to the simulator. These steps are implemented in lines 690-770 of "srs/main.cpp" and conclude the implementation of the path planner.

## Testing

The method described above had been simulated using the Udacity simulator. The vehicle was capable of driving at least 30 Miles without any incidents while maintaining a speed close to the reference velocity of 50 mph as seen in the picture below:

![test_pic][test]

## Frenet Coordinates

Frenet coordinates are heavily used throughout my implementation of the code. The basic concept of frenet coordinates is very simple with the "s" values representing the distance traveled along the road regardless of its shape while the "d" value represents the distance from the center of the road. A simple visualization on the differance between cartesian coordinates (left picture) and frenet coordinates (right picture) can be seen below:

![frenet_cor][frenet]

