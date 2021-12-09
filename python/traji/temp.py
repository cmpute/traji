import _lib2
import numpy as np

path = _lib2.Path(np.array([(0, 1), (2, 3)]))
print(path)
for p in path:
    print(p)
print(np.array(path))
print(path.numpy())
print(path.respacing(0.02))

traj = _lib2.Trajectory(path, [1, 2])
print(traj.point_at(_lib2.PathPosition(0, 0)))
print(len(traj))
print(traj.numpy())

print(traj.velocity_at(_lib2.PathPosition(0, 0)))
