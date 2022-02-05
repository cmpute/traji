import traji
import numpy as np

def test_quintpoly():
    # simple case
    x0 = [0, 1, 0]
    xT = [1, 0, 0]
    y0 = [1, 0, 0]
    yT = [0, 1, 0]
    qtraj = traji.QuinticPolyTrajectory(5, x0, xT, y0, yT)

    assert qtraj.T == 5
    assert len(qtraj.x_coeffs) == 6
    assert len(qtraj.y_coeffs) == 6

    traj = qtraj.periodize(0.1)
    arr = traj.numpy()
    assert np.all(arr[:, 2] >= 0) and np.all(arr[:, 2] <= 5)
    assert len(traj) == len(arr) - 1
    assert np.all(np.diff(arr[:, 2]) <= 0.1 + 1e-5)
    for i, p in enumerate(arr):
        assert np.allclose(qtraj.point_at(p[2]), p[:2])

    # create random quintic polynomial
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
        assert np.allclose(traj.acceleration_at(T), [xT[2], yT[2]], atol=5e-3, rtol=5e-3)

def test_trajectory_interpolation():
    x0 = [0, 1, 0]
    xT = [1, 0, 0]
    y0 = [1, 0, 0]
    yT = [0, 1, 0]
    traj = traji.QuinticPolyTrajectory(5, x0, xT, y0, yT)
    traj = traj.periodize(0.02)

    for s in np.linspace(-0.1,1.1,13) * traj.length:
        traj.point_at(s)
        traj.velocity_at(s)
        traj.acceleration_at(s)

def test_ctratrajectory():
    traj = traji.CTRATrajectory(5, traji.Point(0, 0), 1, 1, 1, 1)
    traj.point_at(1)
    
    assert traj.T == 5
    traj = traj.periodize(0.1)

    assert traj.timestamps[-1] == 5

def test_velocity():
    traj = traji.Trajectory(traji.Path([(0, 0), (3, 4), (6, 8), (9, 12)]), [0, 2, 7, 9])
    speed52 = [3/2, 4/2]
    assert np.allclose(traj.velocity_from(2.5), speed52)
    assert np.allclose(traj.velocity_at(1), speed52)
    assert np.allclose(traj.velocity_from(12.5), speed52)
    assert np.allclose(traj.velocity_at(8), speed52)

    assert np.allclose(traj.velocity_from(2.5, interpolate=True), speed52)
    assert np.allclose(traj.velocity_at(1, interpolate=True), speed52)
    assert np.allclose(traj.velocity_from(12.5, interpolate=True), speed52)
    assert np.allclose(traj.velocity_at(8, interpolate=True), speed52)

if __name__ == "__main__":
    test_quintpoly()
    test_trajectory_interpolation()
    test_ctratrajectory()
    test_velocity()
