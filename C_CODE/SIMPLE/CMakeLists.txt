cmake_minimum_required(VERSION 3.15...4.0)

project(
	SIMPLECavity
	VERSION 0.1
	LANGUAGES CXX)

add_executable(
	SIMPLE
	main.cpp 
	MomentumEquations.cpp 
	Corrections.cpp 
	Boundary.cpp 
	SIMPLE.cpp
)
# target_include_directories(SIMPLE PRIVATE /home/wujorgen/eigen-3.4.0)
if (APPLE)
	include_directories(/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/)
elseif(UNIX)
	include_directories(/usr/local/include/eigen-3.4.0)
endif()

target_compile_options(SIMPLE PRIVATE -march=native)

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
