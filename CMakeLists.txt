# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project(Rod)

cmake_minimum_required(VERSION 2.8.11)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

# CGAL and its components
find_package(CGAL QUIET COMPONENTS)

if (NOT CGAL_FOUND)

    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()

endif ()

# include helper file
include(${CGAL_USE_FILE})


# Boost and its components
find_package(Boost REQUIRED)

if (NOT Boost_FOUND)

    message(STATUS "This project requires the Boost library, and will not be compiled.")

    return()

endif ()

# include for local directory

# include for local package


# Creating entries for target: executable
# ############################

add_executable(QueryHandler QueryHandler.cpp NaiveQueryHandler.cpp
               NaiveQueryHandler.h IQueryHandler.h CGAL_defines.h
               SegmentCheck.h SegmentCheck.cpp DijkstraComputePath.cpp
               MyQueryHandler.cpp MyQueryHandler.h)
add_executable(RodPathFinder RodPathFinder.cpp MyRodPathFinder.cpp MyRodPathFinder.h Path.cpp Path.h NaiveQueryHandler.cpp NaiveQueryHandler.h IQueryHandler.h CGAL_defines.h)
add_to_cached_list(CGAL_EXECUTABLE_TARGETS QueryHandler)
add_to_cached_list(CGAL_EXECUTABLE_TARGETS RodPathFinder)

# Link the executable to CGAL and third-party libraries
target_link_libraries(QueryHandler ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
target_link_libraries(RodPathFinder ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})