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
#include "spline.h"

using namespace std;
using namespace std::chrono;
// driving parameters

// speed limit
const double max_speed_mph = 49.5;

// time in between 2 path points (determined by simulator)
const double time_between_path_points = 0.02;

// maximum acceleration in m/s2
const double max_acceleration_mps2 = 5; // typical acceleration for a car
const double max_deceleration_mps2 = 7; // typical deceleration for a car
const double max_acceleration_pp = max_acceleration_mps2 * time_between_path_points;
const double max_deceleration_pp = max_deceleration_mps2 * time_between_path_points;

// distance to car ahead that triggers a lane change (m)
const double overtake_trigger_dist = 30;

// required free space needed in adjacent lane to conduct lane change (m)
const double overtake_ahead_dist_required = 40; // ahead of car
const double overtake_behind_dist_required = -5; // behind car

// minimum time our car has to stay in a line before it can make a lane change to 5 s
const milliseconds min_time_in_lane_ms = milliseconds(5000);


const double move_right_min_speed = 40;
const double move_right_ahead_dist_required = 60; // ahead of car
const double move_right_behind_dist_required = -30; // behind car

const double max_speed_tracking_distance = 60;


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
  
  double heading = atan2( (map_y-y),(map_x-x) );
  
  double angle = abs(theta-heading);
  
  if(angle > pi()/4)
  {
    closestWaypoint++;
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

/*   Returns whether a car is within a certain distance of current car in a given lane
 *
 *   Parameters:
 *     - sensor_fusion = data structure contain information on the surrounding cars
 *     - prev_size = number of points in the car's path sent to the simulator
 *     - lane = lane in which to look for a car [0, 1, or 2]
 *     - car_s = the car's s Frenet coordinate (distance along the road) [m]
 *     - look_behind, look_ahead = distances behind and ahead of the car to look for a car [m]
 */
bool car_is_in_lane(vector<vector<double>> sensor_fusion, int prev_size, int lane, double car_s,
                    double look_behind, double look_ahead) {
  bool result = false;
  for (int i=0; i < sensor_fusion.size(); i++) {
    float d = sensor_fusion[i][6];
    // the -3/+3 deviate from the project's Q&A. I added another additional 1 m buffer to account for cars
    // that are starting to change lanes
    if (d > (2 + 4 * lane - 3) && d < (2 + 4 * lane + 3)) {
      double vx = sensor_fusion[i][3];
      double vy = sensor_fusion[i][4];
      double check_speed = sqrt(vx * vx + vy * vy);
      double check_car_s = sensor_fusion[i][5];
      
      check_car_s += ((double)prev_size * time_between_path_points * check_speed);
      
      if ((check_car_s - car_s > look_behind) && (check_car_s-car_s < look_ahead)) {
        result = true;
        
      }
    }
  }
  return result;
}

/*   Returns the speed of a car at a certain distance ahead of the car's current position.
 *
 *   Parameters:
 *     - sensor_fusion = data structure contain information on the surrounding cars
 *     - prev_size = number of points in the car's path sent to the simulator
 *     - lane = lane in which to look for a car [0, 1, or 2]
 *     - car_s = the car's s Frenet coordinate (distance along the road) [m]
 *     - look_behind, look_ahead = distances behind and ahead of the car to look for a car [m]
 */
double get_speed_of_car_ahead(vector<vector<double>> sensor_fusion, int prev_size, int lane, double car_s,
                    double look_behind, double look_ahead) {
  double result = -1;
  for (int i=0; i < sensor_fusion.size(); i++) {
    float d = sensor_fusion[i][6]; // distance to the side of the road
    // the -3/+3 deviate from the project's Q&A. I added another additional 1 m buffer to account for cars
    // that are starting to change lanes
    if (d > (2 + 4 * lane - 3) && d < (2 + 4 * lane + 3)) {
      double vx = sensor_fusion[i][3];
      double vy = sensor_fusion[i][4];
      double check_speed = sqrt(vx * vx + vy * vy);
      double check_car_s = sensor_fusion[i][5];
      
      check_car_s += ((double)prev_size * time_between_path_points * check_speed);
      
      if ((check_car_s - car_s > look_behind) && (check_car_s - car_s < look_ahead)) {
        result = check_speed * 2.24; // conversion of m/s to mph
      }
    }
  }
  return result;
}
/*   Calculate the car's speed adjusted to the speed of the car ahead taking into account the
 *   car's maximum acceleration/deceleration parameters.
 *
 *   Parameters:
 *     - current_speed = current car speed [mph]
 *     - ahead_speed = speed of the car ahead [mph]
 */
double adjust_speed_to_car_ahead(double current_speed, double ahead_speed) {
  if (ahead_speed < 0) ahead_speed = max_speed_mph;
  if (ahead_speed > max_speed_mph) ahead_speed = max_speed_mph;
  double target_speed = min(max_speed_mph, ahead_speed);
  double new_speed = 0;
  if (current_speed > ahead_speed) {
    new_speed = max(current_speed - max_deceleration_pp, target_speed);
  } else {
    new_speed = min(current_speed + max_acceleration_pp, target_speed);
  }
  return new_speed;
}

int main() {
  uWS::Hub h;
  
  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;
  
  // Waypoint map to read from
  string map_file_ = "data/highway_map.csv";
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
  
  int lane = 1;
  double reference_velocity = 0.; // mph
  milliseconds last_lane_change_ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
  
  
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &lane, &reference_velocity, &last_lane_change_ms](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
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
          
          int prev_size = previous_path_x.size();
          
          if (prev_size > 0) {
            car_s = end_path_s;
          }
          
          bool too_close = false;
          
          // Is there a car ahead in the lane that is too close?
          too_close = car_is_in_lane(sensor_fusion, prev_size, lane, car_s, 0, overtake_trigger_dist);
          
          milliseconds curr_time_ms = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
          
          // keep track of the time since the last lane change. Only after a set minimum time since the last
          // lane change is exceeded, will we make another lane changes to avoid making too many lan changes
          bool min_time_in_lane_exceeded = (curr_time_ms - last_lane_change_ms) > min_time_in_lane_ms;
          
          // The actual state transition model
          if (too_close) {
            if (lane == 1) { // the car ahead is too close and we're in the middle lane
              // if there is space in the left lane and the minimum time before making another lane change is exceeded,
              // go to the left lane
              
              // note that we are verifying that there is space by looking at the 5m behind our position and
              // ahead_dist_required (40m) in front. The behind_dist_required (-5m) is intended to detect cars that
              // are going faster than us but who are still behind
              if (!car_is_in_lane(sensor_fusion, prev_size, 0, car_s,
                                  overtake_behind_dist_required,
                                  overtake_ahead_dist_required) &&
                  min_time_in_lane_exceeded) {
                lane = 0;
                last_lane_change_ms = curr_time_ms;
              // if there is no space in the left lane, but there is space in the right lane and minimum time is exceeded
              // change lanes to the right
              } else if (!car_is_in_lane(sensor_fusion, prev_size, 2, car_s,
                                         overtake_behind_dist_required,
                                         overtake_ahead_dist_required) &&
                         min_time_in_lane_exceeded) {
                lane = 2;
                last_lane_change_ms = curr_time_ms;
              } else {
                // there is no space in the left or right lanes, so adjust speed to the speed of the car ahead
                double ahead_speed = get_speed_of_car_ahead(sensor_fusion, prev_size, lane, car_s, 0,
                                                            overtake_trigger_dist);
                reference_velocity = adjust_speed_to_car_ahead(reference_velocity, ahead_speed);
              }
            } else if (lane == 2) { // the car ahead is too close and we're in the right lane
              // if there is space in middle lane, and minimum time is exceeded, change lanes to left
              if (!car_is_in_lane(sensor_fusion, prev_size, 1, car_s,
                                  overtake_behind_dist_required,
                                  overtake_ahead_dist_required) &&
                  min_time_in_lane_exceeded) {
                lane = 1;
                last_lane_change_ms = curr_time_ms;
              } else {
                // else, stay in lane and adjust speed to car ahead
                double ahead_speed = get_speed_of_car_ahead(sensor_fusion, prev_size, lane, car_s, 0,
                                                            overtake_trigger_dist);
                reference_velocity = adjust_speed_to_car_ahead(reference_velocity, ahead_speed);
              }
            } else if (lane == 0) { // the car ahead is too close and we're in the left lane
              // if there is space in middle lane, and minimum time is exceeded, change lanes to left
              if (!car_is_in_lane(sensor_fusion, prev_size, 1, car_s,
                                  overtake_behind_dist_required,
                                  overtake_ahead_dist_required) &&
                  min_time_in_lane_exceeded) {
                lane = 1;
                last_lane_change_ms = curr_time_ms;
              } else {
                // else, stay in lane and adjust speed to car ahead
                double ahead_speed = get_speed_of_car_ahead(sensor_fusion, prev_size, lane, car_s, 0,
                                                            overtake_trigger_dist);
                reference_velocity = adjust_speed_to_car_ahead(reference_velocity, ahead_speed);
              }
            }
          } else { // normal driving conditions, we're not too close to car ahead
            // if there is a car ahead a bit further, get its speed to make sure we're not accelerating
            // too much beyond that car's speed to avoid having to break excessively in front of the car
            double ahead_speed = get_speed_of_car_ahead(sensor_fusion, prev_size, lane, car_s, 0,
                                                        max_speed_tracking_distance);
            // allow our car to drive a little faster than a car 60m ahead, so we can eventually take over
            if (ahead_speed > 0) ahead_speed = ahead_speed + 5;
            reference_velocity = adjust_speed_to_car_ahead(reference_velocity, ahead_speed);
            // if there is enough space, no car coming, and we have a minimum speed, change to the rightmost lane
            // but one lane at a time
            if (lane < 2 && !car_is_in_lane(sensor_fusion, prev_size, lane + 1, car_s,
                                            move_right_behind_dist_required,
                                            move_right_ahead_dist_required) &&
                min_time_in_lane_exceeded && reference_velocity > move_right_min_speed) {
              lane += 1;
              last_lane_change_ms = curr_time_ms;
            }
          }
          
          // ptsx and ptsy are vectors of cartesian coordinates that we will use to create the path out of through
          // spline interpolation
          vector<double> spline_pts_x;
          vector<double> spline_pts_y;
          
          double ref_x, ref_y, ref_yaw;
          double prev_x, prev_y;
          
          // First, add the reference and previous positions to spline_pts_x and spline_pts_y
          // The reference position (ref_x, ref_y) is the position from where we need to calculate new path points
          // The previous point (prev_x, prev_y) is the point before the reference position.
          if (prev_size < 2) {
            // In the beginning (when no previous path is returned from the simulator (i.e., prev_size < 2),
            // the reference position is the current position of the car and the previous position is calculated
            // based on the car's direction/ yaw.
            ref_x = car_x;
            ref_y = car_y;
            ref_yaw = deg2rad(car_yaw);
            // The 5 is a semi-arbitrary distance we assume for the previous position. The exact value
            // does not matter, as long as the direction is correct
            prev_x = ref_x - 5 * cos(ref_yaw);
            prev_y = ref_y - 5 * sin(ref_yaw);
          } else {
            // After 1 timestep, the simulator returns a previous path with processed points removed, and we
            // need to only add 1-2 new points at the back of the simulator.
            // Hence, the reference position is the last element of the previous path
            ref_x = previous_path_x[prev_size - 1];
            ref_y = previous_path_y[prev_size - 1];
            prev_x = previous_path_x[prev_size - 2];
            prev_y = previous_path_y[prev_size - 2];
            ref_yaw = atan2(ref_y - prev_y, ref_x - prev_x);
          }
          
          spline_pts_x.push_back(prev_x);
          spline_pts_x.push_back(ref_x);
          spline_pts_y.push_back(prev_y);
          spline_pts_y.push_back(ref_y);
          
          // To calculate the rest of the path, get the cartesian coordinates of 3 waypoints that are
          // 30, 60, and 90m ahead of the car's current position but in the same lane (Frenet coordinates)
          // This also assumes that a lane change takes 30m, which is reasonable.
          // if the first distance is set too low (i.e., 10m), then the lane change occurs too abrupt
          // if it is set too high, lane changes occur too slowly
          vector<double> next_wp0 = getXY(car_s + 30, 2 + 4 * lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s + 60, 2 + 4 * lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s + 90, 2 + 4 * lane, map_waypoints_s, map_waypoints_x, map_waypoints_y);
          
          spline_pts_x.push_back(next_wp0[0]);
          spline_pts_x.push_back(next_wp1[0]);
          spline_pts_x.push_back(next_wp2[0]);
          
          spline_pts_y.push_back(next_wp0[1]);
          spline_pts_y.push_back(next_wp1[1]);
          spline_pts_y.push_back(next_wp2[1]);
          
          // transform coordinates from cartesian to car XY coordinates in-place
          for (int i=0; i < spline_pts_x.size(); i++) {
            double shift_x = spline_pts_x[i] - ref_x;
            double shift_y = spline_pts_y[i] - ref_y;
            
            spline_pts_x[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
            spline_pts_y[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
          }
          
          // create a spline based on these waypoints
          tk::spline s;
          s.set_points(spline_pts_x, spline_pts_y);
          
          vector<double> new_path_x;
          vector<double> new_path_y;
          
          for (int i = 0; i < prev_size; i++) {
            new_path_x.push_back(previous_path_x[i]);
            new_path_y.push_back(previous_path_y[i]);
          }
          
          // define a path made up of (x,y) points that the car will visit sequentially every 20 ms
          for (int i = 1; i <= 50 - prev_size; i++) {
            // car moves to a new waypoint every 20ms, calculate the distance a car would drive in 20 ms
            double distance = (0.02 * reference_velocity / 2.24);
            // For each additional path point needed, calculate car XY position if car would drive this distance
            // along the X axis (i.e., along its current yaw) for 1 unit of 20 ms
            double pt_car_coord_x = i * distance;
            double pt_car_coord_y = s(pt_car_coord_x);
            
            // transform from car XY coordinates to map XY coordinates
            double pt_map_coord_x = pt_car_coord_x * cos(ref_yaw) - pt_car_coord_y * sin(ref_yaw);
            double pt_map_coord_y = pt_car_coord_x * sin(ref_yaw) + pt_car_coord_y * cos(ref_yaw);
            pt_map_coord_x += ref_x;
            pt_map_coord_y += ref_y;
            
            new_path_x.push_back(pt_map_coord_x);
            new_path_y.push_back(pt_map_coord_y);
          }
          
          json msgJson;
          msgJson["next_x"] = new_path_x;
          msgJson["next_y"] = new_path_y;
          
          auto msg = "42[\"control\","+ msgJson.dump()+"]";
          
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
