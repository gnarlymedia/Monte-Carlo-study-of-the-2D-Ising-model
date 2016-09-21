cmake_minimum_required(VERSION 3.6)

project(proj_4)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(SOURCE_FILES_1 proj_4.c)
add_executable(proj_4 ${SOURCE_FILES_1})

#project(tute9)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set(SOURCE_FILES_2 tute9.c)
#add_executable(tute9 ${SOURCE_FILES_2})

find_library(M_LIB m X11)
find_library(CPGPLOT_LIB cpgplot pgplot png z)
target_link_libraries(proj_4 ${M_LIB})
target_link_libraries(proj_4 ${CPGPLOT_LIB})