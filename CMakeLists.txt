cmake_minimum_required(VERSION 2.8.3)

project(utils)

# Add support for c++11
add_definitions(-std=c++11 -O3)

# Find required packages
find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy)


find_package(Eigen 3 QUIET) 
if(NOT Eigen_FOUND)
	include(cmake/FindEigen3.cmake)
endif() 

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS roscpp rospy
 )

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(include ${catkin_INCLUDE_DIRS})

FILE(GLOB SRCS src/Utils.cpp)
FILE(GLOB INCS include/Utils.h)


add_library(${PROJECT_NAME} ${SRCS} ${INCS})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})


add_executable(test_utils src/test_utils.cpp)
target_link_libraries(test_utils ${catkin_LIBRARIES} ${PROJECT_NAME})
add_dependencies(test_utils ${PROJECT_NAME}_gencpp)