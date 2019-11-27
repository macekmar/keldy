#  Copyright Olivier Parcollet 2014
#  Copyright Simons Foundation 2019
#    Author: Nils Wentzell, Philipp Dumitrescu

#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

#
# This module looks for cuba.
# It sets up : CUBA_INCLUDE_DIR, CUBA_LIBRARIES
# Use CUBA_ROOT to specify a particular location
#

if(CUBA_INCLUDE_DIR AND CUBA_LIBRARIES)
  set(CUBA_FIND_QUIETLY TRUE)
endif()

find_path(CUBA_INCLUDE_DIR
  NAMES cuba.h
  PATHS
    ${CUBA_ROOT}/include
    $ENV{CUBA_ROOT}/include
    ENV CPATH
    ENV C_INCLUDE_PATH
    ENV CPLUS_INCLUDE_PATH
    ENV OBJC_INCLUDE_PATH
    ENV OBJCPLUS_INCLUDE_PATH
    /usr/include
    /usr/local/include
    /opt/local/include
    /sw/include
  DOC "Include Directory for Cuba"
)

find_library(CUBA_LIBRARIES
  NAMES cuba
  PATHS
    ${CUBA_INCLUDE_DIR}/../lib
    ${CUBA_ROOT}/lib
    $ENV{CUBA_ROOT}/lib
    ENV LIBRARY_PATH
    ENV LD_LIBRARY_PATH
    /usr/lib
    /usr/local/lib
    /opt/local/lib
    /sw/lib
  DOC "Cuba library"
)

include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CUBA DEFAULT_MSG CUBA_LIBRARIES CUBA_INCLUDE_DIR)

mark_as_advanced(CUBA_INCLUDE_DIR CUBA_LIBRARIES)

# file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/cuba_api_test.cpp "#include \"cuba.h\"\nint main(){ cuba_plan p; return p.cuba_flags; }")
# try_compile(CUBA_OLD_API ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/cuba_api_test.cpp COMPILE_DEFINITIONS -I${CUBA_INCLUDE_DIR} )
# file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/a.out ${CMAKE_CURRENT_BINARY_DIR}/cuba_api_test.cpp)

# Interface target
# We refrain from creating an imported target since those cannot be exported
add_library(cuba INTERFACE)
target_link_libraries(cuba INTERFACE ${CUBA_LIBRARIES})
target_include_directories(cuba SYSTEM INTERFACE ${CUBA_INCLUDE_DIR})
