CC=/illumina/thirdparty/gcc/gcc-4.9.2/bin/g++ -O3 -std=c++11 -static

BOOST_BASE=/illumina/thirdparty/boost/boost_1_54_0_python2.7/

BOOST_INC_DIR=${BOOST_BASE}/include
BOOST_LIB_DIR=${BOOST_BASE}/lib
BOOST_LIB_SPEC=-lboost_program_options -lboost_filesystem -lboost_system -lboost_regex -lboost_date_time -lpthread

SAMTOOLS_BASE=~/popdata/progs/htslib/
SAMTOOLS_INC_DIR=${SAMTOOLS_BASE}
SAMTOOLS_LIB_DIR=${SAMTOOLS_BASE}
# actually -lhts -lz BUT -lz already specified separately
SAMTOOLS_LIB_SPEC=-lhts

COMPR_LIB_SPEC=-lz

SRC_DIR=src
INC_DIR=-I. -Ithird_party
LIB_DIR=src
INSTALL_DIR=bin
BUILD_DIR=local-build

LIBS=${BUILD_DIR}/region.o ${BUILD_DIR}/parameters.o \
     ${BUILD_DIR}/purity.o ${BUILD_DIR}/pair_collector.o \
     ${BUILD_DIR}/ref_genome.o ${BUILD_DIR}/bam_file.o \
     ${BUILD_DIR}/irr_finder.o ${BUILD_DIR}/common.o \
     ${BUILD_DIR}/classifiers.o

all: ${INSTALL_DIR} ${INSTALL_DIR}/ExpansionHunterDenovo

${BUILD_DIR}:
	mkdir -p ${BUILD_DIR}

${BUILD_DIR}/%.o: common/%.cc
	${CC} -c -o $@ $< ${INC_DIR}

${BUILD_DIR}/%.o: purity/%.cc
	${CC} -c -o $@ $< ${INC_DIR}

${BUILD_DIR}/%.o: regions/%.cc
	${CC} -c -o $@ $< -I ${BOOST_INC_DIR} ${INC_DIR}

${BUILD_DIR}/%.o: classification/%.cc
	${CC} -c -o $@ $< -I ${BOOST_INC_DIR} ${INC_DIR}

${BUILD_DIR}/%.o: pair_collector/%.cc
	${CC} -c -o $@ $< ${INC_DIR}

${BUILD_DIR}/%.o: ${LIB_DIR}/%.cc ${BUILD_DIR}
	${CC} -c -o $@ $< -I ${SAMTOOLS_INC_DIR} -I ${BOOST_INC_DIR} ${INC_DIR}

${INSTALL_DIR}:
	mkdir -p ${INSTALL_DIR}

${INSTALL_DIR}/ExpansionHunterDenovo: ${SRC_DIR}/expansion_hunter_denovo.cc ${LIBS}
	${CC} -o $@ $^ -I ${SAMTOOLS_INC_DIR} -I ${BOOST_INC_DIR} ${INC_DIR} -L ${SAMTOOLS_LIB_DIR} -L ${BOOST_LIB_DIR} ${SAMTOOLS_LIB_SPEC} ${COMPR_LIB_SPEC} ${BOOST_LIB_SPEC}
