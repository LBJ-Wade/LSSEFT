cmake_minimum_required(VERSION 3.2)
project(LSSEFT)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Ofast")

SET ( SPLINTER_LIBRARY "splinter-static-2-0" )

find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS system filesystem timer date_time mpi program_options serialization)

# pull in bundled SPLINTER sources; on most platforms will need to switch default build to 64-bit
set ( ARCH "x86-64" )
ADD_SUBDIRECTORY ( "thirdparty/splinter" "splinter-build" )

# use version of Eigen3 bundled with SPLINTER
SET ( EIGEN3_INCLUDE_DIR  "thirdparty/splinter/thirdparty/Eigen" )

include_directories(./
  thirdparty/splinter/include
  ${EIGEN3_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${MPI_CXX_INCLUDE_PATH})

set(LOCALIZATION
  localizations/messages.h)

link_directories("/Users/ds283/Documents/Code/Cuba/Cuba-4.2")
link_directories("/Users/david/Documents/Code/Cuba/Cuba-4.2")

set(LOCALIZATION_EN_GB
  localizations/en_GB/en_GB.h
  localizations/en_GB/command_line.h
  localizations/en_GB/error_handler.h
  localizations/en_GB/task_manager.h
  localizations/en_GB/database.h
  localizations/en_GB/sqlite3.h
  localizations/en_GB/transaction_manager.h
  localizations/en_GB/range.h
  localizations/en_GB/format.h
  localizations/en_GB/master_controller.h
  localizations/en_GB/tree_power_spectrum.h)

set(TOP_LEVEL_SOURCE_FILES
  main.cpp
  core.h
  defaults.h
  exceptions.h)

set(CONTROLLER_SOURCE_FILES
  controller/argument_cache.cpp controller/argument_cache.h
  controller/task_manager.cpp controller/task_manager.h
  controller/master_controller.cpp controller/master_controller.h
  controller/slave_controller.cpp controller/slave_controller.h
  controller/local_environment.cpp controller/local_environment.h
  controller/scheduler.cpp controller/scheduler.h)

set(UTILITY_SOURCE_FILES
  utilities/python_finder.cpp utilities/python_finder.h
  utilities/formatter.cpp utilities/formatter.h)

set(ERROR_SOURCE_FILES
  error/error_handler.cpp error/error_handler.h)

set(MISC_SOURCE_FILES
  misc/ansi_colour_codes.h)

set(MPI_SOURCE_FILES
  MPI_detail/mpi_operations.h)

set(COSMOLOGY_SOURCE_FILES
  cosmology/FRW_model.cpp cosmology/FRW_model.h
  cosmology/Planck_defaults.h
  cosmology/types.h
  cosmology/concepts/range.h
  cosmology/concepts/range_detail/stepping.h
  cosmology/concepts/range_detail/range.h
  cosmology/concepts/range_detail/common.h
  cosmology/concepts/range_detail/aggregation.h
  cosmology/transfer_integrator.cpp cosmology/transfer_integrator.h
  cosmology/concepts/oneloop.cpp cosmology/concepts/oneloop.h
  cosmology/concepts/transfer_function.cpp cosmology/concepts/transfer_function.h
  cosmology/oneloop_integrator.cpp cosmology/oneloop_integrator.h
  cosmology/constants.h
  cosmology/concepts/tree_power_spectrum.cpp cosmology/concepts/tree_power_spectrum.h)

set(UNITS_SOURCE_FILES
  units/eV_units.h)

set(DATABASE_SOURCE_FILES
  database/data_manager.cpp database/data_manager.h
  database/tokens.cpp database/tokens.h
  database/transaction_manager.cpp database/transaction_manager.h
  database/wavenumber_database.h
  database/z_database.cpp database/z_database.h
  database/generic_record_iterator.h database/generic_value_iterator.h
  database/powerspectrum_database.cpp database/powerspectrum_database.h
  database/z_record.cpp database/z_record.h
  database/wavenumber_record.h
  database/powerspectrum_record.cpp database/powerspectrum_record.h
  database/k_database.h database/IR_database.h database/UV_database.h)

set(SQLITE3_DETAIL_SOURCE_FILES
  sqlite3_detail/operations.h
  sqlite3_detail/sqlite3_defaults.h
  sqlite3_detail/create.cpp sqlite3_detail/create.h
  sqlite3_detail/utilities.cpp sqlite3_detail/utilities.h
  sqlite3_detail/FRW_model.cpp sqlite3_detail/FRW_model.h
  sqlite3_detail/sqlite3_policy.cpp sqlite3_detail/sqlite3_policy.h
  sqlite3_detail/redshift.cpp sqlite3_detail/redshift.h
  sqlite3_detail/wavenumber.h
  sqlite3_detail/temporary_tables.cpp sqlite3_detail/temporary_tables.h
  sqlite3_detail/missing_elements.cpp sqlite3_detail/missing_elements.h
  sqlite3_detail/store.cpp sqlite3_detail/store.h)

set(SOURCE_FILES
  ${TOP_LEVEL_SOURCE_FILES}
  ${CONTROLLER_SOURCE_FILES}
  ${UTILITY_SOURCE_FILES}
  ${ERROR_SOURCE_FILES}
  ${MPI_SOURCE_FILES}
  ${MISC_SOURCE_FILES}
  ${COSMOLOGY_SOURCE_FILES}
  ${UNITS_SOURCE_FILES}
  ${DATABASE_SOURCE_FILES}
  ${SQLITE3_DETAIL_SOURCE_FILES}
  ${LOCALIZATION}
  ${LOCALIZATION_EN_GB})

add_executable(lsseft
  main.cpp
  controller/argument_cache.cpp
  controller/local_environment.cpp
  controller/master_controller.cpp
  controller/slave_controller.cpp
  controller/task_manager.cpp
  controller/scheduler.cpp
  error/error_handler.cpp
  utilities/python_finder.cpp
  utilities/formatter.cpp
  cosmology/FRW_model.cpp
  cosmology/concepts/transfer_function.cpp
  cosmology/concepts/oneloop.cpp
  cosmology/concepts/tree_power_spectrum.cpp
  cosmology/transfer_integrator.cpp
  cosmology/oneloop_integrator.cpp
  database/data_manager.cpp
  database/tokens.cpp
  database/transaction_manager.cpp
  database/z_database.cpp
  database/powerspectrum_database.cpp
  database/z_record.cpp
  database/powerspectrum_record.cpp
  sqlite3_detail/utilities.cpp
  sqlite3_detail/create.cpp
  sqlite3_detail/FRW_model.cpp
  sqlite3_detail/redshift.cpp
  sqlite3_detail/temporary_tables.cpp
  sqlite3_detail/missing_elements.cpp
  sqlite3_detail/store.cpp
  sqlite3_detail/sqlite3_policy.cpp)

target_link_libraries(lsseft sqlite3 cuba ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${SPLINTER_LIBRARY})

add_executable(dummy_clion_target EXCLUDE_FROM_ALL ${SOURCE_FILES})
