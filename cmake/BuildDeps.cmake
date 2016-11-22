CMAKE_MINIMUM_REQUIRED(VERSION 3.0)

INCLUDE(ExternalProject)

FUNCTION(BUILD_SPLINTER BUILD_DIR INSTALL_PREFIX)

  SET(SPLINTER_SRC ${BUILD_DIR}/src)
  SET(SPLINTER_BIN ${BUILD_DIR}/bin)
  SET(SPLINTER_INSTALL ${INSTALL_PREFIX})
  SET(SPLINTER_INSTALL ${INSTALL_PREFIX} PARENT_SCOPE)

  EXTERNALPROJECT_ADD(
    SPLINTER
    GIT_REPOSITORY https://github.com/bgrimstad/splinter.git
    GIT_TAG 2c877c0d7f68025764f7af048cc609e5c05080e0
    SOURCE_DIR ${SPLINTER_SRC}
    BINARY_DIR ${SPLINTER_BIN}
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=${SPLINTER_INSTALL} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  )

  SET(SPLINTER_LIBRARIES ${SPLINTER_INSTALL}/lib/${CMAKE_FIND_LIBRARY_PREFIXES}splinter-static-3-0${CMAKE_STATIC_LIBRARY_SUFFIX} PARENT_SCOPE)
  SET(SPLINTER_INCLUDE_DIRS ${SPLINTER_INSTALL}/include ${SPLINTER_INSTALL}/include/SPLINTER PARENT_SCOPE)

  ADD_DEPENDENCIES(DEPS SPLINTER)

ENDFUNCTION(BUILD_SPLINTER)
