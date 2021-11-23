First build the traji C++ library
```
mkdir build
cd build
cmake ..
cmake --build . --config Release -- -j $(nproc)
```
and then execute `python setup.py install`.
