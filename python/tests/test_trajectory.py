import traji
import numpy as np

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
    test_quintpoly()
