mkdir build
cd build
cmake -DCMAKE_C_COMPILER=//opt/homebrew/Cellar/gcc/14.2.0/bin/gcc-14 -DCMAKE_CXX_COMPILER=/opt/homebrew/Cellar/gcc/14.2.0/bin/g++-14 -DCMAKE_BUILD_TYPE=Release ..
cmake  --build . --parallel
cmake --install .
cd ..