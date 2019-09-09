add_library(mergeworkflow STATIC
        MergeWorkflow.hh MergeWorkflow.cpp
        MergeParameters.hh MergeParameters.cpp
        MultisampleProfile.hh MultisampleProfile.cpp)

target_link_libraries(mergeworkflow io Boost::filesystem region)
target_include_directories(mergeworkflow PUBLIC ${CMAKE_SOURCE_DIR})
