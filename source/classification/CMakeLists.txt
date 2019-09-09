add_library(classification irr_finder.cc classifiers.cc)
target_link_libraries(classification purity common)
add_subdirectory(unit_tests)