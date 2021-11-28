import traji

def test_project():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    print(path.project(traji.Point(0, 0.1)))
    print(path.project(traji.Point(0.1, 0)))

def test_respacing():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    print(path.respacing(0.1))
    print(path.respacing(0.1, smooth_radius=0.3))
    print(path.respacing(0.3))
    print(path.respacing(0.3, smooth_radius=0.3))

def test_densify():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    print(path.densify(0.1))
    print(path.densify(0.3))

def test_smooth():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    print(path.smooth(0.1))
    print(path.smooth(0.1, smooth_radius=0.3))
    print(path.smooth(0.3))
    print(path.smooth(0.3, smooth_radius=0.3))

if __name__ == "__main__":
    test_project()
    test_respacing()
    test_densify()
    test_smooth()
