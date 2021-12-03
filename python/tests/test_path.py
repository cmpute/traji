import traji
import numpy as np
from shapely.geometry import Point as sPoint, LineString

def test_project():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    d, pos = path.project(traji.Point(0, 0.1))
    assert np.isclose(d, 0)
    assert pos.segment == 0 and np.isclose(pos.fraction, 0.1)
    d, pos = path.project(traji.Point(0.1, 0))
    assert np.isclose(d, -0.1)
    assert pos.segment == 0 and np.isclose(pos.fraction, 0)

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

def test_shapely_interop():
    point = traji.Point([1, 2])
    shapely_point = sPoint(point)
    shapely_point = point.shapely()
    point2 = traji.Point(shapely_point)
    assert point == point2

    path = traji.Path([(0, 0), (1, 2), (2, 1)])
    shapely_path = path.shapely()
    path2 = traji.Path(shapely_path)
    assert path == path2

if __name__ == "__main__":
    # test_project()
    # test_respacing()
    # test_densify()
    test_shapely_interop()
