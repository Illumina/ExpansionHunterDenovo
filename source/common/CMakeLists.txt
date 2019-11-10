add_library(common STATIC
        Parameters.hh
        SequenceUtils.hh SequenceUtils.cpp Interval.cpp Interval.hh)

target_include_directories(common PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(common Boost::boost)

#add_executable(SequenceUtilsTest SequenceUtilsTest.cpp)
#target_link_libraries(SequenceUtilsTest common)
#add_test(NAME SequenceUtilsTest COMMAND SequenceUtilsTest)
