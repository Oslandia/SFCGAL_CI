zypper update -y
zypper install -y \
    cmake cgal-devel ninja gcc gcc-c++ \
    libboost_serialization-devel \
    libboost_program_options-devel \
    libboost_test-devel \
    nlohmann_json-devel
rm -rf build
