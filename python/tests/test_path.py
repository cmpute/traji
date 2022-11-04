from curses import window
import traji
import numpy as np
from shapely.geometry import Point as sPoint, LineString

def test_deduplicate():
    path = traji.Path([(0, 0), (0, 0), (1, 1), (1, 1)])
    assert path == traji.Path([(0, 0), (1, 1)])

def test_project():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    d, pos = path.project(traji.Point(0, 0.1))
    assert np.isclose(d, 0)
    assert pos.segment == 0 and np.isclose(pos.fraction, 0.1)
    d, pos = path.project(traji.Point(0.1, 0))
    assert np.isclose(d, -0.1)
    assert pos.segment == 0 and np.isclose(pos.fraction, 0)

def test_position_conversion():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    for s in np.linspace(-0.1,1.1,13) * path.length:
        pos = traji.PathPosition.from_s(path, s)

        if s >= 0 and s <= path.length:
            assert 0 <= pos.fraction and pos.fraction < 1
        assert np.isclose(pos.to_s(path), s)

def test_respacing():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    respaced = path.respacing(0.1)
    assert np.allclose(respaced.segment_lengths, np.linalg.norm(np.diff(respaced, axis=0), axis=1))
    assert np.all(np.array(respaced.segment_lengths) < 0.1 + 1e-5)
    respaced = path.respacing(0.1, smooth_radius=0.3)
    assert np.allclose(respaced.segment_lengths, np.linalg.norm(np.diff(respaced, axis=0), axis=1))
    assert np.all(np.array(respaced.segment_lengths) < 0.1 + 1e-5)
    respaced = path.respacing(0.3)
    assert np.allclose(respaced.segment_lengths, np.linalg.norm(np.diff(respaced, axis=0), axis=1))
    assert np.all(np.array(respaced.segment_lengths) < 0.3 + 1e-5)
    respaced = path.respacing(0.3, smooth_radius=0.3)
    assert np.allclose(respaced.segment_lengths, np.linalg.norm(np.diff(respaced, axis=0), axis=1))
    assert np.all(np.array(respaced.segment_lengths) < 0.3 + 1e-5)

def test_densify():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    densified = path.densify(0.1)
    assert np.all(np.array(densified.segment_lengths) < 0.1 + 1e-5)
    densified = path.densify(0.3)
    assert np.all(np.array(densified.segment_lengths) < 0.3 + 1e-5)

def test_resample():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])

    s_list = np.random.rand(10) * 3 - 0.5
    assert path.resample_from(s_list) == traji.Path([path.point_from(s) for s in s_list])

    s_list = np.sort(np.random.rand(10) * 3 - 0.5)
    assert path.resample_from(s_list) == traji.Path([path.point_from(s) for s in s_list])

def test_shapely_interop():
    point = traji.Point([1, 2])
    shapely_point = sPoint(point.x, point.y)
    shapely_point = point.shapely()
    point2 = traji.Point(shapely_point)
    assert point == point2

    path = traji.Path([(0, 0), (1, 2), (2, 1)])
    shapely_path = path.shapely()
    path2 = traji.Path(shapely_path)
    assert path == path2

def test_curvature():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])

    assert path.curvature_from(-0.5) == 0
    assert path.curvature_from(0) == 0
    assert np.isclose(path.curvature_from(1), -np.sqrt(2))
    assert path.curvature_from(2) == 0
    assert path.curvature_from(2.5) == 0

    # draw a circle
    grid = np.linspace(-np.pi, np.pi, 100)
    points = np.array([np.cos(grid), np.sin(grid)]).T
    path = traji.Path(points)
    np.isclose(path.curvature_from(1), 1)
    np.isclose(path.curvature_from(2), 1)

def test_extract():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])

    assert path.extract_from(0.5, 1.5) == traji.Path([(0, 0.5), (0, 1), (0.5, 1)])
    assert path.extract_from(-2, -1) == traji.Path([(0, -2), (0, -1)])
    assert path.extract_from(3, 4) == traji.Path([(2, 1), (3, 1)])
    assert path.extract_from(-2, 4) == traji.Path([(0, -2), (0, 1), (3, 1)])
    assert path.extract_from(-1, 0.5) == traji.Path([(0, -1), (0, 0.5)])
    assert path.extract_from(1.5, 3) == traji.Path([(0.5, 1), (2, 1)])

if __name__ == "__main__":
    test_project()
    test_respacing()
    test_densify()
    test_resample()
    test_shapely_interop()
    test_position_conversion()
    test_curvature()
    test_extract()
