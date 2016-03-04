cmake_minimum_required(VERSION 2.8)
project(correctReadcount)

SET(CMAKE_BUILD_TYPE Release)

add_subdirectory(lib)

include_directories (
	${correctReadcount_SOURCE_DIR}/lib/split
	${correctReadcount_SOURCE_DIR}/lib/fastahack
	${correctReadcount_SOURCE_DIR}/lib/kent/inc
	${correctReadcount_SOURCE_DIR}/lib/bamtools
)

add_executable (gcCounter src/bin/gcCounter.cpp)
add_executable (mapCounter src/bin/mapCounter.cpp)
add_executable (readCounter src/bin/readCounter.cpp)

target_link_libraries(gcCounter fastahack)
target_link_libraries(mapCounter kent split)
target_link_libraries(readCounter bamtools split)

SET_TARGET_PROPERTIES(gcCounter PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/bin)
SET_TARGET_PROPERTIES(mapCounter PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/bin)
SET_TARGET_PROPERTIES(readCounter PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${correctReadcount_BINARY_DIR}/bin)

add_subdirectory(src/util)
