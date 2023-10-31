sudo yum update -qy
sudo yum install -y \
    cmake boost boost-devel gmp gmp-c++ gmp-devel mpfr mpfr-devel make \
    xz

#CGAL
wget https://github.com/CGAL/cgal/releases/download/v"$1"/CGAL-"$1".tar.xz
tar xJf CGAL-"$1".tar.xz
cd CGAL-"$1" && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$CI_PROJECT_DIR/CGAL .. && make && make install && cd ../..

cmake --version
