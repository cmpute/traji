First build the traji C++ library
```shell
mkdir build
cd build
cmake ..
cmake --build . --config Release -- -j $(nproc)
```
and then execute `python setup.py install`.
