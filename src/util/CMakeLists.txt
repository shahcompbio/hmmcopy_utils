include_directories (
	${correctReadcount_SOURCE_DIR}/lib/split
	${correctReadcount_SOURCE_DIR}/lib/fastahack
	${correctReadcount_SOURCE_DIR}/lib/kent/inc
)

add_executable (bigWigInfo bigwig/bigWigInfo.c)
add_executable (bigWigSummary bigwig/bigWigSummary.c)
add_executable (bigWigToBedGraph bigwig/bigWigToBedGraph.c)
add_executable (bigWigToWig bigwig/bigWigToWig.c)
add_executable (wigToBigWig bigwig/wigToBigWig.c)

add_executable (segToMap seg/segToMap.cpp)
add_executable (segToGc seg/segToGc.cpp)

add_executable (fastaToRead mappability/internal/fastaToRead.cpp)

target_link_libraries(bigWigInfo kent)
target_link_libraries(bigWigSummary kent)
target_link_libraries(bigWigToBedGraph kent)
target_link_libraries(bigWigToWig kent)
target_link_libraries(wigToBigWig kent)

target_link_libraries(segToMap kent split)
target_link_libraries(segToGc fastahack)

target_link_libraries(fastaToRead fastahack)

SET_TARGET_PROPERTIES(bigWigInfo PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/util/bigwig)
SET_TARGET_PROPERTIES(bigWigSummary PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/util/bigwig)
SET_TARGET_PROPERTIES(bigWigToBedGraph PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/util/bigwig)
SET_TARGET_PROPERTIES(bigWigToWig PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/util/bigwig)
SET_TARGET_PROPERTIES(wigToBigWig PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/util/bigwig)

SET_TARGET_PROPERTIES(segToMap PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/util/seg)
SET_TARGET_PROPERTIES(segToGc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/util/seg)

SET_TARGET_PROPERTIES(fastaToRead PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/util/mappability/internal)

file(COPY ${correctReadcount_SOURCE_DIR}/src/util/mappability/internal/readToMap.pl DESTINATION ${correctReadcount_BINARY_DIR}/util/mappability/internal)
file(COPY ${correctReadcount_SOURCE_DIR}/src/util/mappability/generateMap.pl DESTINATION ${correctReadcount_BINARY_DIR}/util/mappability)
file(COPY ${correctReadcount_SOURCE_DIR}/src/util/renameChr.pl DESTINATION ${correctReadcount_BINARY_DIR}/util/)
