add_library(outlierworkflow STATIC
        OutlierWorkflow.hh OutlierWorkflow.cpp
        OutlierParameters.hh OutlierParameters.cpp
        ZScoreAnalysis.hh ZScoreAnalysis.cpp)

target_link_libraries(outlierworkflow io Boost::filesystem region common)
target_include_directories(outlierworkflow PUBLIC ${CMAKE_SOURCE_DIR})
