cd /home/gmehdevi/sgoinfre
wget https://www.vtk.org/files/release/9.3/VTK-9.3.0.tar.gz
tar -xzvf VTK-9.3.0.tar.gz
mkdir VTK-9.3.0-build
cd VTK-9.3.0-build
cmake ../VTK-9.3.0 -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)