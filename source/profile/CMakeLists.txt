add_library(profileworkflow STATIC
        ProfileWorkflow.hh ProfileWorkflow.cpp
        ProfileParameters.hh ProfileParameters.cpp
        SampleRunStats.hh SampleRunStats.cpp)

target_link_libraries(profileworkflow io Boost::filesystem)
target_include_directories(profileworkflow PUBLIC ${CMAKE_SOURCE_DIR})
