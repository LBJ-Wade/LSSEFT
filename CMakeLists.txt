cmake_minimum_required(VERSION 3.2)
project(LSSEFT)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem timer date_time mpi program_options)

include_directories(${Boost_INCLUDE_DIRS}
                    ${MPI_CXX_INCLUDE_PATH})

set(SOURCE_FILES main.cpp)

add_executable(lsseft ${SOURCE_FILES})

target_link_libraries(lsseft ${Boost_LIBRARIES})
