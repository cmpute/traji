# TrajI
A C++/Python library designed for trajectory calculation in Cartesian / Frenet coordinate system

## Basic Usage
This library include several basic concepts implemented as classes:
- `traji.Point`: A 2D point in Cartesian or Frenet coordinate
- `traji.Path`: A line string consists of a list of consecutive points.
- `traji.Trajectory`: A line string path with associated timestamps for each point.
- `traji.HeteroPath`: A path consists of heterogenerous segments, including line segments, arc segments and other curve types.
- `traji.PathPosition`: This class define a position along a path. It can be converted to a distance value or a timestamp.
- `traji.QuinticPolyTrajectory`: An optimal trajectory defined as a quintic polynomial with regard to time.
- `traji.CTRATrajectory`: A trajectory generated from constant turn-rate and acceleration movement.

The library contains the following functionalities for these concepts with efficient implementations:
- Transformation between the concepts (e.g. `Path` from/to `HeteroPath`, `Path` from/to `Trajectory`)
- Transformation between Cartesian coordinate and Frenet coordinate
- Indexing on the linear concepts (`Point` from/to `PathPosition`)
- Binary operations between the concepts (e.g. intersection, distance)

For C++ API, please refer to the single header file [traji.hpp](include/traji.hpp) for the library.

For Python API, please refer to the type annotation file [_lib.pyi](python/traji/_lib.pyi) for an overview of the definitions.

## Reference

- Frenet coordinates: https://github.com/fjp/frenet
- MATLAB traj library: https://www.mathworks.com/help/nav/ref/trajectoryoptimalfrenet.cart2frenet.html
- Bezier library: https://github.com/oysteinmyrmo/bezier
