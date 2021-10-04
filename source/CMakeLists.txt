cmake_minimum_required(VERSION 3.10)
project(ExpansionHunterDenovo LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

enable_testing()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

add_subdirectory(thirdparty/boost-cmake)

# Required to make spdlog work
include_directories(thirdparty)

add_subdirectory(common)
add_subdirectory(io)
add_subdirectory(reads)
add_subdirectory(region)
add_subdirectory(profile)
add_subdirectory(merge)

add_executable(ExpansionHunterDenovo
        app/ExpansionHunterDenovo.cpp
        app/GatherParameters.hh
        app/GatherParameters.cpp
        app/Version.hh)
target_include_directories(ExpansionHunterDenovo PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(ExpansionHunterDenovo PRIVATE
        Boost::system Boost::program_options
        common region profileworkflow mergeworkflow)

add_executable(UnitTests
        tests/UnitTests.cpp
        tests/SequenceUtilsTest.cpp
        tests/PurityScoreTest.cpp
        tests/GenomicRegionTest.cpp
        tests/IrrFinderTest.cpp)
target_link_libraries(UnitTests common reads region)
target_include_directories(UnitTests PUBLIC ${CMAKE_SOURCE_DIR})


#list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
#
#enable_testing()
#
#include(ExternalProject)
#
########################## Google Test ############################
## Download and unpack googletest at configure time
#configure_file(cmake/google_test.cmake googletest-download/CMakeLists.txt)
#execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
#		RESULT_VARIABLE result
#		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
#if(result)
#	message(FATAL_ERROR "CMake step for googletest failed: ${result}")
#endif()
#execute_process(COMMAND ${CMAKE_COMMAND} --build .
#		RESULT_VARIABLE result
#		WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download )
#if(result)
#	message(FATAL_ERROR "Build step for googletest failed: ${result}")
#endif()
#
## Add googletest directly to our build. This defines
## the gtest and gtest_main targets.
#add_subdirectory(${CMAKE_BINARY_DIR}/googletest-main
#		${CMAKE_BINARY_DIR}/googletest-build)
###################################################################
#
#ExternalProject_Add(zlib
#	PREFIX ${CMAKE_BINARY_DIR}/thirdparty/zlib
#	GIT_REPOSITORY "https://github.com/madler/zlib.git"
#	GIT_TAG "v1.2.8"
#	UPDATE_COMMAND ""
#	BUILD_IN_SOURCE 1
#	CONFIGURE_COMMAND ${CMAKE_BINARY_DIR}/thirdparty/zlib/main/zlib/configure --prefix=${CMAKE_BINARY_DIR}/thirdparty/zlib --static
#	INSTALL_DIR ${CMAKE_BINARY_DIR}/thirdparty/zlib
#	LOG_DOWNLOAD 1
#	LOG_INSTALL 1
#)
#
#ExternalProject_Add(htslib
#	PREFIX ${CMAKE_BINARY_DIR}/thirdparty/htslib
#        GIT_REPOSITORY "https://github.com/samtools/htslib.git"
#	GIT_TAG "1.3.1"
#        UPDATE_COMMAND ""
#	BUILD_IN_SOURCE 1
#        CONFIGURE_COMMAND ""
#	BUILD_COMMAND make
#	INSTALL_COMMAND make install prefix=${CMAKE_BINARY_DIR}/thirdparty/htslib
#	LOG_DOWNLOAD 1
#)
#
#include_directories(${CMAKE_BINARY_DIR}/thirdparty/zlib/include)
#set(zlib_static ${CMAKE_BINARY_DIR}/thirdparty/zlib/lib/libz.a)
#set(htslib_static ${CMAKE_BINARY_DIR}/thirdparty/htslib/lib/libhts.a)
#
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#
#find_package(Boost 1.4 COMPONENTS program_options filesystem regex date_time system REQUIRED)
#
#include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories(${Boost_INCLUDE_DIR})
#include_directories(${CMAKE_BINARY_DIR}/thirdparty/htslib/include)
#include_directories(${CMAKE_SOURCE_DIR}/thirdparty)
#
#add_subdirectory(common)
#add_subdirectory(purity)
#add_subdirectory(regions)
#add_subdirectory(classification)
#add_subdirectory(reads)
#
#file(GLOB SOURCES "main/*.cc")

#add_executable(ExpansionHunterDenovo ${SOURCES})
#target_compile_features(ExpansionHunterDenovo PRIVATE cxx_range_for)

#add_dependencies(htslib zlib)
#add_dependencies(ExpansionHunterDenovo htslib)

#target_link_libraries(ExpansionHunterDenovo regions purity classification reads pthread ${htslib_static} ${zlib_static} ${Boost_LIBRARIES})
