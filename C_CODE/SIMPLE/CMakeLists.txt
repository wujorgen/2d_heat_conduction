cmake_minimum_required(VERSION 3.15...4.0)

project(
	SIMPLECavity
	VERSION 0.1
	LANGUAGES CXX)

add_executable(SIMPLECavity main.cpp momentum_eqn.cpp corrections.cpp boundary.cpp)
# target_include_directories(SIMPLECavity PRIVATE /home/wujorgen/eigen-3.4.0)
if (APPLE)
	include_directories(/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/)
elseif(UNIX)
	include_directories(/usr/local/include/eigen-3.4.0)
endif()

target_compile_options(SIMPLECavity PRIVATE -march=native)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
