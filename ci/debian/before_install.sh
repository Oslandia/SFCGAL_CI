export DEBIAN_FRONTEND=noninteractive
apt-get update -qq
apt-get install --yes \
    cmake libboost-chrono-dev libboost-program-options-dev libboost-timer-dev \
    libboost-test-dev libboost-thread-dev nlohmann-json3-dev \
    libboost-system-dev libboost-serialization-dev \
    libmpfr-dev libgmp-dev \
    xz-utils
#CGAL

wget https://github.com/CGAL/cgal/releases/download/v"$1"/CGAL-"$1".tar.xz
tar xJf CGAL-"$1".tar.xz
cd CGAL-"$1" && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$CI_PROJECT_DIR/CGAL .. && make && make install && cd ../..

cmake --version
