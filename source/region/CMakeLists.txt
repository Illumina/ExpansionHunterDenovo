add_library(region STATIC
        GenomicRegion.hh GenomicRegion.cpp
        ReferenceContigInfo.hh ReferenceContigInfo.cpp)
target_include_directories(region PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(region Boost::boost)
