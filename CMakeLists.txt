CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
PROJECT(LSSEFT)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/)
INCLUDE(BuildDeps)

# find MPI libraries
FIND_PACKAGE(MPI REQUIRED)

# find Boost libraries
FIND_PACKAGE(Boost REQUIRED COMPONENTS system filesystem timer date_time mpi program_options serialization)

# find OpenSSL libraries
FIND_PACKAGE(OpenSSL REQUIRED)

# pull in bundled SPLINTER sources; on most platforms will need to switch default build to 64-bit
SET(ARCH "x86-64")

ADD_CUSTOM_TARGET(DEPS ALL)
SET(EXTERNAL_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/thirdparty/dist")

BUILD_SPLINTER(${CMAKE_BINARY_DIR}/thirdparty/SPLINTER ${EXTERNAL_INSTALL_PREFIX})

# use version of Eigen3 bundled with SPLINTER
SET(EIGEN3_INCLUDE_DIR "thirdparty/splinter/thirdparty/Eigen")

# Cuba libraries and include files
SET(CUBA_INCLUDE_DIR "/Users/ds283/Documents/Code/Cuba/Cuba-4.2" "/Users/david/Documents/Code/Cuba/Cuba-4.2" "/usr/local/include")
SET(CUBA_LIBRARY_DIR "/Users/ds283/Documents/Code/Cuba/Cuba-4.2" "/Users/david/Documents/Code/Cuba/Cuba-4.2" "/usr/local/lib")

SET(LOCALIZATION
  localizations/messages.h)

LINK_DIRECTORIES(${CUBA_LIBRARY_DIR})

SET(LOCALIZATION_EN_GB
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

SET(TOP_LEVEL_SOURCE_FILES
  main.cpp
  core.h
  defaults.h
  exceptions.h)

SET(CONTROLLER_SOURCE_FILES
  controller/argument_cache.cpp controller/argument_cache.h
  controller/task_manager.cpp controller/task_manager.h
  controller/master_controller.cpp controller/master_controller.h
  controller/slave_controller.cpp controller/slave_controller.h
  controller/local_environment.cpp controller/local_environment.h
  controller/scheduler.cpp controller/scheduler.h)

SET(UTILITY_SOURCE_FILES
  utilities/python_finder.cpp utilities/python_finder.h
  utilities/formatter.cpp utilities/formatter.h)

SET(ERROR_SOURCE_FILES
  error/error_handler.cpp error/error_handler.h)

SET(MISC_SOURCE_FILES
  misc/ansi_colour_codes.h)

SET(MPI_SOURCE_FILES
  MPI_detail/mpi_operations.h
  MPI_detail/mpi_traits.h
  MPI_detail/mpi_payloads.cpp MPI_detail/mpi_payloads.h)

SET(COSMOLOGY_SOURCE_FILES
  cosmology/FRW_model.cpp cosmology/FRW_model.h
  cosmology/Planck_defaults.h
  cosmology/types.h
  cosmology/constants.h
  cosmology/concepts/range.h
  cosmology/concepts/range_detail/stepping.h
  cosmology/concepts/range_detail/range.h
  cosmology/concepts/range_detail/common.h
  cosmology/concepts/range_detail/aggregation.h
  cosmology/concepts/oneloop_growth.cpp cosmology/concepts/oneloop_growth.h
  cosmology/concepts/transfer_function.cpp cosmology/concepts/transfer_function.h
  cosmology/concepts/power_spectrum_detail/generic.h
  cosmology/concepts/power_spectrum_detail/types.h
  cosmology/concepts/power_spectrum_detail/linear.cpp cosmology/concepts/power_spectrum_detail/linear.h
  cosmology/concepts/power_spectrum.h
  cosmology/concepts/multipole_Pk.cpp cosmology/concepts/multipole_Pk.h
  cosmology/concepts/oneloop_Pk.cpp cosmology/concepts/oneloop_Pk.h
  cosmology/concepts/Matsubara_A.cpp cosmology/concepts/Matsubara_A.h
  cosmology/oneloop_integrands/shared.h
  cosmology/oneloop_integrands/delta_13.h
  cosmology/oneloop_integrands/delta_22.h
  cosmology/oneloop_integrands/delta_rsd_13.h
  cosmology/oneloop_integrands/delta_rsd_22.h
  cosmology/oneloop_momentum_integrator.cpp cosmology/oneloop_momentum_integrator.h
  cosmology/concepts/loop_integral.cpp cosmology/concepts/loop_integral.h
  cosmology/transfer_integrator.cpp cosmology/transfer_integrator.h
  cosmology/oneloop_growth_integrator.cpp cosmology/oneloop_growth_integrator.h
  cosmology/oneloop_Pk_calculator.cpp cosmology/oneloop_Pk_calculator.h
  cosmology/multipole_Pk_calculator.cpp cosmology/multipole_Pk_calculator.h)

SET(UNITS_SOURCE_FILES
  units/Mpc_units.h)

SET(DATABASE_SOURCE_FILES
  database/data_manager.cpp database/data_manager.h database/data_manager_impl.h
  database/tokens.cpp database/tokens.h
  database/transaction_manager.cpp database/transaction_manager.h
  database/wavenumber_database.h
  database/z_database.cpp database/z_database.h
  database/generic_record_iterator.h database/generic_value_iterator.h
  database/powerspectrum_database.h
  database/z_record.cpp database/z_record.h
  database/wavenumber_record.h
  database/powerspectrum_record.h
  database/k_database.h
  database/IR_cutoff_database.h
  database/UV_cutoff_database.h
  database/IR_resum_database.h)

SET(SQLITE3_DETAIL_SOURCE_FILES
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
  sqlite3_detail/store.cpp sqlite3_detail/store.h
  sqlite3_detail/find.cpp sqlite3_detail/find.h
  sqlite3_detail/power_spectrum.cpp sqlite3_detail/power_spectrum.h)

SET(SOURCE_FILES
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

ADD_EXECUTABLE(lsseft
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
  cosmology/concepts/oneloop_growth.cpp
  cosmology/concepts/loop_integral.cpp
  cosmology/concepts/oneloop_Pk.cpp
  cosmology/concepts/multipole_Pk.cpp
  cosmology/concepts/Matsubara_A.cpp
  cosmology/concepts/power_spectrum_detail/linear.cpp
  cosmology/transfer_integrator.cpp
  cosmology/oneloop_growth_integrator.cpp
  cosmology/oneloop_momentum_integrator.cpp
  cosmology/oneloop_Pk_calculator.cpp
  cosmology/multipole_Pk_calculator.cpp
  database/data_manager.cpp
  database/tokens.cpp
  database/transaction_manager.cpp
  database/z_database.cpp
  database/z_record.cpp
  MPI_detail/mpi_payloads.cpp
  sqlite3_detail/utilities.cpp
  sqlite3_detail/create.cpp
  sqlite3_detail/FRW_model.cpp
  sqlite3_detail/redshift.cpp
  sqlite3_detail/temporary_tables.cpp
  sqlite3_detail/missing_elements.cpp
  sqlite3_detail/store.cpp
  sqlite3_detail/find.cpp
  sqlite3_detail/power_spectrum.cpp
  sqlite3_detail/sqlite3_policy.cpp)

ADD_DEPENDENCIES(lsseft DEPS)

TARGET_LINK_LIBRARIES(lsseft sqlite3 cuba ${Boost_LIBRARIES} ${MPI_LIBRARIES} ${SPLINTER_LIBRARIES} ${OPENSSL_LIBRARIES})

TARGET_COMPILE_OPTIONS(lsseft PRIVATE -std=c++14)

TARGET_INCLUDE_DIRECTORIES(lsseft PRIVATE
  ./
  ${CMAKE_CURRENT_BINARY_DIR}
  ${SPLINTER_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${OPENSSL_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_PATH}
  ${CUBA_INCLUDE_DIR})

ADD_EXECUTABLE(dummy_clion_target EXCLUDE_FROM_ALL ${SOURCE_FILES})
