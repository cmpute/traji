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

    for _ in range(50):
        x0 = np.random.rand(3)
        xT = np.random.rand(3)
        y0 = np.random.rand(3)
        yT = np.random.rand(3)
        T = np.random.rand() * 10
        traj = traji.QuinticPolyTrajectory(T, x0, xT, y0, yT)

        assert np.allclose(traj.point_at(0).numpy(), [x0[0], y0[0]])
        assert np.allclose(traj.point_at(T).numpy(), [xT[0], yT[0]], atol=1e-3, rtol=1e-3)
        assert np.allclose(traj.velocity_at(0), [x0[1], y0[1]])
        assert np.allclose(traj.velocity_at(T), [xT[1], yT[1]], atol=1e-3, rtol=1e-3)
        assert np.allclose(traj.acceleration_at(0), [x0[2], y0[2]])
        assert np.allclose(traj.acceleration_at(T), [xT[2], yT[2]], atol=1e-3, rtol=1e-3)

if __name__ == "__main__":
    test_quintpoly()
