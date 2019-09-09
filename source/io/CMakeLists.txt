find_package(Threads)

include(ExternalProject)
include(FindZLIB)
include(FindBZip2)
include(FindLibLZMA)

#ExternalProject_Add(htslib
#        PREFIX ${CMAKE_BINARY_DIR}/thirdparty/htslib
#        GIT_REPOSITORY "https://github.com/samtools/htslib.git"
#        GIT_TAG "1.3.1"
#        UPDATE_COMMAND ""
#        BUILD_IN_SOURCE 1
#        CONFIGURE_COMMAND "./configure"
#        BUILD_COMMAND make
#        INSTALL_COMMAND make install prefix=${CMAKE_BINARY_DIR}/thirdparty/htslib
#        LOG_DOWNLOAD 1
#        )


ExternalProject_Add(htslib_project
	URL ${CMAKE_SOURCE_DIR}/thirdparty/htslib-1.9.tar.bz2
	SOURCE_DIR ${CMAKE_BINARY_DIR}/thirdparty/htslib-1.9-source
	INSTALL_DIR ${CMAKE_BINARY_DIR}/thirdparty/htslib-1.9
	CONFIGURE_COMMAND ""
	BUILD_COMMAND make -C <SOURCE_DIR> prefix=<INSTALL_DIR>
	INSTALL_COMMAND make -C <SOURCE_DIR> install prefix=<INSTALL_DIR>
	)

#set(htslib_static ${CMAKE_BINARY_DIR}/thirdparty/htslib/lib/libhts.a)
add_library(htslib STATIC IMPORTED)
add_dependencies(htslib htslib_project)
set_property(TARGET htslib 
	PROPERTY IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/thirdparty/htslib-1.9/lib/libhts.a)


add_library(io STATIC
        HtsFileStreamer.hh HtsFileStreamer.cpp
        HtsHelpers.hh HtsHelpers.cpp
        Reference.hh Reference.cpp)

target_include_directories(io PUBLIC
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_BINARY_DIR}/thirdparty/htslib-1.9/include
        )

target_link_libraries(io PUBLIC
        reads
        htslib
 	${ZLIB_LIBRARIES}
        ${LIBLZMA_LIBRARIES}
	${BZIP2_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
	)

add_dependencies(io htslib)

