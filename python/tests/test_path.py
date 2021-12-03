import traji
import numpy as np

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

def test_quintpoly():
    x0 = [0, 1, 0]
    xT = [1, 0, 0]
    y0 = [1, 0, 0]
    yT = [0, 1, 0]
    traj = traji.QuinticPolyTrajectory(5, x0, xT, y0, yT)

    assert traj.T == 5
    assert len(traj.x_coeffs) == 6
    assert len(traj.y_coeffs) == 6

    traj = traj.periodize(0.1)
    arr = traj.numpy()
    assert np.all(arr[:, 2] >= 0) and np.all(arr[:, 2] <= 5)
    assert len(traj) == len(arr)
    assert np.all(np.diff(arr[:, 2]) <= 0.1 + 1e-5)

if __name__ == "__main__":
    test_project()
    test_respacing()
    test_densify()
    test_quintpoly()
