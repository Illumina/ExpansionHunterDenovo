add_library(reads STATIC
        Read.hh Read.cpp
        ../profile/PairCollector.hh ../profile/PairCollector.cpp
        ../profile/ReadClassification.hh ../profile/ReadClassification.cpp
        IrrFinder.hh IrrFinder.cpp
        Purity.hh Purity.cpp)
target_link_libraries(reads common)
