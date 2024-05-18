cmake_minimum_required(VERSION 3.27)
cmake_policy(SET CMP0057 NEW)
set(CMAKE_CXX_FLAGS "-std=c++11")

project(BijectiveDiscretizedRotations)

#Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

add_definitions(-g)

aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} ${DGTAL_LIBRARIES})
