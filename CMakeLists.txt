cmake_minimum_required(VERSION 3.2)
project(LSSEFT)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem timer date_time mpi program_options)

include_directories(./
                    ${Boost_INCLUDE_DIRS}
                    ${MPI_CXX_INCLUDE_PATH})

set(SOURCE_FILES main.cpp controller/argument_cache.cpp controller/task_manager.cpp controller/master_controller.cpp controller/slave_controller.cpp controller/local_environment.cpp utilities/python_finder.cpp error/error_handler.cpp)

add_executable(lsseft ${SOURCE_FILES})

target_link_libraries(lsseft ${Boost_LIBRARIES} ${MPI_LIBRARIES})
