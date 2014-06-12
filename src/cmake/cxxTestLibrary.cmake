#
# Manta
# Copyright (c) 2013-2014 Illumina, Inc.
#
# This software is provided under the terms and conditions of the
# Illumina Open Source Software License 1.
#
# You should have received a copy of the Illumina Open Source
# Software License 1 along with this program. If not, see
# <https://github.com/sequencing/licenses/>
#

################################################################################
##
## Configuration file for the unit tests subdirectory
##
## author Ole Schulz-Trieglaff
##
################################################################################

set(IS_QUIET true)
include(${THIS_CXX_EXECUTABLE_CMAKE})

file (GLOB TEST_SOURCE "*.cpp")

set(TEST_TARGET unit_test_${THIS_LIB_DIR})

add_executable(${TEST_TARGET} ${TEST_SOURCE})
add_dependencies(${TEST_TARGET} STARKA_OPT)

if (THIS_LIBRARY_SOURCES)
    set(ADDITIONAL_UNITTEST_LIB ${ADDITIONAL_UNITTEST_LIB} starka_${THIS_LIB_DIR})
endif ()

target_link_libraries (${TEST_TARGET} ${ADDITIONAL_UNITTEST_LIB} ${THIS_AVAILABLE_LIBRARIES}
                      ${SAMTOOLS_LIBRARY} ${TABIX_LIBRARY} ${Boost_LIBRARIES} ${THIS_ADDITIONAL_LIB})

set(TEST_BINARY ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARGET})

add_test(${TEST_TARGET} ${TEST_BINARY} "--log_level=test_suite")
add_dependencies(STARKA_UNITTESTS ${TEST_TARGET})