cmake_minimum_required(VERSION 3.15...4.0)

project(
	ChorinProjectionCavity
	VERSION 0.1
	LANGUAGES CXX)

add_executable(ChorinProjectionCavity main.cpp)
# target_include_directories(ChorinProjectionCavity PRIVATE /home/wujorgen/eigen-3.4.0)
if (APPLE)
	include_directories(/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/)
elseif(UNIX)
	include_directories(/usr/local/include/eigen-3.4.0)
endif()

target_compile_options(ChorinProjectionCavity PRIVATE -march=native)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
