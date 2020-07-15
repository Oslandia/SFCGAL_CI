export DEBIAN_FRONTEND=noninteractive
sudo apt-get update -qq
sudo apt-get install --yes \
    cmake libboost-chrono-dev libboost-program-options-dev libboost-filesystem-dev libboost-timer-dev \
    libboost-test-dev libboost-thread-dev \
    libboost-system-dev libboost-serialization-dev \
    libmpfr-dev libgmp-dev \
    cmake
#CGAL

wget https://github.com/CGAL/cgal/releases/download/releases/CGAL-"$1"/CGAL-"$1".tar.xz
tar xJf CGAL-"$1".tar.xz
cd CGAL-"$1" && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$CI_PROJECT_DIR/CGAL .. && make && make install && cd ../..

cmake --version
