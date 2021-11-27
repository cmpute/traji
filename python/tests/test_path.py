import traji

def test_project():
    path = traji.Path([(0, 0), (0, 1), (1, 1)])
    print(path.project(traji.Point(0, 0.1)))
    print(path.project(traji.Point(0.1, 0)))

if __name__ == "__main__":
    test_project()