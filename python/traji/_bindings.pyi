from enum import Enum
from typing import Union, Tuple, Optional, List

class SegmentType(Enum):
    Line
    Arc
    QuadraticBezier
    CubicBezier
    Polynomial

class Point:
    def __init__(self,
                 x: Union[float, Tuple[float, float]],
                 y: Optional[float] = None) -> None: ...
    @property
    def x(self) -> float: ...
    @property
    def y(self) -> float: ...
    def __array__(self) -> 'numpy.ndarray': ...
    def numpy(self) -> 'numpy.ndarray': ...
    def shapely(self) -> 'shapely.geometry.Point': ...

class PathPosition:
    def __init__(self,
                 component: Optional[int] = 0,
                 segment: Optional[int] = 0,
                 fraction: Optional[float] = 0): ...
    @property
    def component(self) -> int: ...
    @property
    def segment(self) -> int: ...
    @property
    def fraction(self) -> float: ...
    @staticmethod
    def from_s(path: Path, s: Union[float, List[float]]) -> Union[PathPosition, List[PathPosition]]: ...
    def to_s(self, path: Path) -> float: ...
    @staticmethod
    def from_t(traj: Trajectory, t: Union[float, List[float]]) -> Union[PathPosition, List[PathPosition]]: ...
    def to_t(self, traj: Trajectory) -> float: ...

class Path:
    def __init__(self, points: List[Union[Point, Tuple[float, float]]]) -> None: ...
    @property
    def length(self) -> float: ...
    @property
    def segment_lengths(self) -> List[float]: ...

    def numpy(self) -> 'numpy.ndarray': ...
    def shapely(self) -> 'shapely.geometry.LineString': ...
    def point_from(self, s: float) -> Point: ...
    def point_at(self, pos: PathPosition) -> Point: ...
    def tangent_from(self, s: float) -> float: ...
    def tangent_at(self, pos: PathPosition) -> float: ...
    def interpolate_from(self, values: List[float], s: float) -> List[float]: ...
    def interpolate_at(self, values: List[float], pos: PathPosition) -> List[float]: ...
    def project(self, point: Point) -> Tuple[float, PathPosition]: ...
    def respacing(self, resolution: float, smooth_radius: float = 0) -> Path: ...
    def densify(self, resolution: float) -> Path: ...
    def smooth(self, smooth_radius: float, segtype: str = "arc") -> HeteroPath: ...
    def resample_from(self, s_list: List[float]) -> Path: ...

class Trajectory(Path):
    def __init__(self,
                 points: Union[Path, List[Union[Point, Tuple[float, float, float]]]],
                 timestamps: Optional[List[float]] = None) -> None: ...
    @property
    def timestamps(self) -> List[float]: ...

    def point_at(self, pos: Union[float, PathPosition]) -> Point: ...
    def tangent_at(self, pos: Union[float, PathPosition]) -> float: ...
    def velocity_from(self, s: float, interpolate: bool = False) -> Tuple[float, float]: ...
    def velocity_at(self, pos: Union[float, PathPosition], interpolate: bool = False) -> Tuple[float, float]: ...
    def acceleration_from(self, s: float, interpolate: bool = False) -> Tuple[float, float]: ...
    def acceleration_at(self, pos: Union[float, PathPosition], interpolate: bool = False) -> Tuple[float, float]: ...
    def resample_at(self, t_list: List[float]): ...

class QuinticPolyTrajectory:
    def __init__(self,
                 T: float,
                 x0: Optional[Tuple[float, float, float]] = None,
                 xT: Optional[Tuple[float, float, float]] = None,
                 y0: Optional[Tuple[float, float, float]] = None,
                 yT: Optional[Tuple[float, float, float]] = None): ...

    @property
    def T(self): ...
    @property
    def x_coeffs(self): ...
    @property
    def y_coeffs(self): ...

    def point_at(self, t: float): ...
    def tangent_at(self, t: float): ...
    def velocity_at(self, t: float): ...
    def acceleration_at(self, t: float): ...
    def periodize(self, interval: float): ...

class CTRATrajectory:
    def __init__(self,
                 T: float, p: Union[List[float], Point],
                 theta: Optional[float] = None,
                 v: Optional[float] = None,
                 a: Optional[float] = None,
                 omega: Optional[float] = None) -> None: ...
    def point_at(self, t: float) -> Point: ...
    def tangent_at(self, t: float) -> float: ...
    def velocity_at(self, t: float) -> float: ...
    def periodize(self, interval: float) -> float: ...

class HeteroPath:
    pass

def frenet_from_cartesian(ref: Path, target: Union[Point, Path, Trajectory]) -> Union[Point, Path, Trajectory]: ...
def frenet_to_cartesian(ref: Path, target: Union[Point, Path, Trajectory]) -> Union[Point, Path, Trajectory]: ...

def distance(lhs: Union[Point, Path], rhs: Union[Point, Path]) -> float: ...
def arg_distance(lhs: Path, rhs: Path) -> Tuple[PathPosition, PathPosition]: ...
def tdistance(lhs: Trajectory, rhs: Trajectory) -> float: ...
def intersection(lhs: Path, rhs: Path) -> List[Point]: ...
def arg_intersection(lhs: Path, rhs: Path) -> List[Tuple[PathPosition, PathPosition]]: ...