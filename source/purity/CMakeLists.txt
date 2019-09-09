add_library(purity purity.cc)
target_link_libraries(purity common)
add_subdirectory(unit_tests)