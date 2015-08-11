cmake_minimum_required(VERSION 3.2)
project(LSSEFT)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem timer date_time mpi program_options)

include_directories(./
  ${Boost_INCLUDE_DIRS}
  ${MPI_CXX_INCLUDE_PATH})

set(LOCALIZATION_EN_GB
  localizations/en_GB/command_line.h
  localizations/en_GB/error_handler.h
  localizations/en_GB/task_manager.h
  localizations/en_GB/en_GB.h)

set(TOP_LEVEL_SOURCE_FILES
  main.cpp
  core.h
  defaults.h)

set(CONTROLLER_SOURCE_FILES
  controller/argument_cache.cpp controller/argument_cache.h
  controller/task_manager.cpp controller/task_manager.h
  controller/master_controller.cpp controller/master_controller.h
  controller/slave_controller.cpp controller/slave_controller.h
  controller/local_environment.cpp controller/local_environment.h)

set(UTILITY_SOURCE_FILES
  utilities/python_finder.cpp utilities/python_finder.h)

set(ERROR_SOURCE_FILES
  error/error_handler.cpp error/error_handler.h)

set(MISC_SOURCE_FILES
  misc/ansi_colour_codes.h)

set(MPI_SOURCE_FILES
  MPI/mpi_operations.h)

set(SOURCE_FILES
  ${TOP_LEVEL_SOURCE_FILES}
  ${CONTROLLER_SOURCE_FILES}
  ${UTILITY_SOURCE_FILES}
  ${ERROR_SOURCE_FILES}
  ${MPI_SOURCE_FILES}
  ${MISC_SOURCE_FILES}
  ${LOCALIZATION_EN_GB})

add_executable(lsseft
  main.cpp
  controller/argument_cache.cpp
  controller/local_environment.cpp
  controller/master_controller.cpp
  controller/slave_controller.cpp
  controller/task_manager.cpp
  error/error_handler.cpp
  utilities/python_finder.cpp)

target_link_libraries(lsseft ${Boost_LIBRARIES} ${MPI_LIBRARIES})

add_executable(dummy_clion_target EXCLUDE_FROM_ALL ${SOURCE_FILES})
