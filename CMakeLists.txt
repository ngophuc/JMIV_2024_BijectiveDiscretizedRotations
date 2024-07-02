cmake_minimum_required(VERSION 3.27)
cmake_policy(SET CMP0057 NEW)
set(CMAKE_CXX_FLAGS "-std=c++11")

project(BijectiveDiscretizedRotations)

#Required in DGtal
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

add_definitions(-g)

#aux_source_directory(. SRC_LIST)
#add_executable(${PROJECT_NAME} ${SRC_LIST})
#target_link_libraries(${PROJECT_NAME} ${DGTAL_LIBRARIES})

add_executable(AlgoBottomUp AlgoBottomUp.cpp CLI11.hpp CircularDoublyLinkedList.h CircularDoublyLinkedList.cpp Tree.h Tree.cpp UtilityFunctions.h UtilityFunctions.cpp Algo1Functions.h Algo1Functions.cpp)
target_link_libraries(AlgoBottomUp ${DGTAL_LIBRARIES})

#add_executable(AlgoTopDown AlgoTopDown.cpp CLI11.hpp DoublyLinkedList.h DoublyLinkedList.cpp Tree.h Tree.cpp)
#target_link_libraries(AlgoTopDown ${DGTAL_LIBRARIES})