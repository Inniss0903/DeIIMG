cmake_minimum_required(VERSION 3.15)
project(DeIIMG)

find_package(OpenCV REQUIRED)
find_package(PkgConfig REQUIRED)

pkg_check_modules(GMP REQUIRED IMPORTED_TARGET gmpxx)

set(CMAKE_CXX_STANDARD 14)

file(GLOB SRC
        "twodimse/*.h"
        "twodimse/*.cpp"
        )

include_directories(twodimse)
include_directories(${OpenCV_INCLUDE_DIRS} PkgConfig::GMP)
add_executable(DeIIMG main.cpp ${SRC} twodimse/imgseg.cpp)

target_link_libraries(DeIIMG PkgConfig::GMP)
target_link_libraries(DeIIMG ${OpenCV_LIBRARIES})
#target_link_directories(DeIIMG PRIVATE /usr/local/lib)
