# Introduction

The model controls the car's speed and yaw through providing a set of (x,y) points to the controller ("the path"). Furthermore, the car assumes a perfect controller and will visit every (x,y) point in the path every 20 ms. Hence, the spacing of the points determines the speed of the car. Following the suggestion of the project walkthrough, we are keeping 50 points in the path, or 1 s of driving. Given a speed of 50 mph, this corresponds to 22 m of distance.

The path planning approach consists of two parts: (1) logic to determine the next points in the path and (2) logic to determine lane changes and speed. We will go through these two aspects individually

# Logic to determine next path points
The controller returns unused path points from the previous iteration. As suggested during the project's Q&A, we add new points to the unprocessed points that were returned from the controller. This ensures smooth transitions.

The new points are calculated by first setting up a spline that follows the direction of the road. The spline waypoints include the car's reference position (either the car's position at the start of the simulation or the last point in the path when the simulation is on-going), its position 20 ms before and 3 points along the road at 30, 60, and 90 m away. This follows the project's Q&A discussion. The position of the first waypoint (30 m) also determines the distance over which a lane change occurs. 30 m seems to be a happy medium to avoid jerk when lane changes are too fast, or swerving behavior when lane changes occur too slow. We use the `getXY` function to transfer coordinates from Frenet (i.e., along the road) to map XY coordinates. Subsequently, the spline waypoint coordinates are transformed from map XY to car XY coordinates to ensure the spline interpolation is well behaved (i.e., only 1 Y coordinate for every X coordinate).

These 5 spline waypoints (i.e., previous position, current position, and positions at 30, 60, and 90 m ahead) are then interpolated to find the additional path points that are missing. At the start of the run, the full path (i.e., 50 points) must be calculated. Once the process has started, we can recycle most of the previously calculated path points and only 1 to 2 points additional points are needed.

Since the car moves to a new waypoint every 20ms, the car's speed determines the distance in between points. Specifically, the distance between points is `0.02 * reference_velocity / 2.24` where reference_velocity is expressed in mph and the 2.24 constant is needed to convert the reference velocity in m/s. Since we are working in car coordinates, this distance is (almost) entirely driven in the X direction, except for minor deviations when the road is curving. Therefore, the XY coordinate of the next path point is simply `(x_ref + distance, spline(x_ref + distance))`, where `spline(...)` is the spline interpolation function. The error that we are making when the road is curving is negligible in practice, as can be observed in the simulator.

Once we have the new waypoints, they are changed back to map XY coordinates and sent to the simulator.

# Logic to determine lane changes and speed

Our logic to determine lane changes and car speed is based on a simple state transition model where the lane is basically the state of the system. While the model is simple and can certainly be improved, I found that it worked well for the conditions of the simulator. It has the following features.
* A car stays in the current lane, unless there is a car less than 30 m ahead that is driving slower than the speed limit (see line 368), there is space in an adjacent lane (i.e., at least 40 m of available space ahead of the car and 5 m behind), and at least 5 s have passed since the last lane change. (logic implemented on line 385 for lane 1, line 407 for lane 2, and line 421 for lane 0)
* Preference is given to take over on the right, but taking over on the left is permitted as well if no space is available to take over on the right.
* If sufficient space is available, the car will try to move to the rightmost (slow) lane, but it will do that 1 lane change at a time. To avoid a situation where we move to the right and then immediately have to take over a car again, the conditions to trigger such a lane change a bit more conservative than for taking over: at least 60 m of space ahead of the car must be available and 30 m of space behind the car in the lane to the right. This logic is implemented on line 444 and further.

In terms of determining the car's speed, we implemented the following logic:
* the default speed is the speed limit `max_speed_mph`. This constant is used in the `adjust_speed_to_car_ahead` function on line 226.
* increasing speed can only occur at a rate determined by the speed's maximum acceleration, which is set to 5 $m/{s^2}$, which is typical for cars. The rate of decreasing speed is set to 7 $m/{s^2}$, which corresponds to going from a speed of 50 mph to standstill in about 3 s. This occurs in the function `adjust_speed_to_car_ahead` on line 226.
* If there is a car less than 60 m ahead in the same lane, then the car will adjust it's speed to no more than 5 mph higher than that car's speed. This adjustment is intended to avoid abrupt speed changes when approaching a slower driving car. This is done on line 437.
* If there is a car less than 30 m ahead which is driving slower and there is no space for a lane change, the vehicles speed will be adjusted to the speed of the car ahead. This is done once for each lane in lines 403, 417, and 431.

# Potential room for improvement

* Appending new points only to unused path-points is not always a good idea since the car will react with a delay of a little under a second to unforeseen circumstances. However, we can probably get away with this since we are driving on a highway, where there is usually good visibility and high predictability of other cars' behavior. As an improvement, the model could determine when a dangerous traffic situation emerges, in which case it would generate a completely new path instantly.
* Instead of triggering lange changes when a large area of available space is available in the adjacent lane, we could have implemented a full trajectory generation approach where we predict the state of all cars around us into the future and only make lane changes when enough space is available given the future state of surrounding cars. This would definitely be safer in more unusual traffic conditions where our car's speed is very different from the speed of the cars around us. This is usually not the case in highway traffic.
* We could have implemented a state transition model with finer-grained states like "preparing for lane change", "tracking speed of car ahead", etc. However, we did not find this necessary to make the car go around the track in the simulator without violating the project's rubric conditions.