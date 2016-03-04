cmake_minimum_required(VERSION 2.8)

# Build the kent library
FILE(GLOB kent_src kent/lib/*.c)
include_directories(kent/inc)
add_library(kent ${kent_src})
target_link_libraries(kent m z)

# Build the split library
include_directories(split)
add_library(split split/split.cpp)

# Build the fastahack library
add_library(fastahack fastahack/Fasta.cpp)
target_link_libraries(fastahack split)

# Build the Bamtools library
FILE(GLOB bamtools_api bamtools/api/*.cpp)
FILE(GLOB bamtools_internal bamtools/api/internal/*.cpp)

# create main BamTools API static library
include_directories(bamtools)
add_library(bamtools ${bamtools_api} ${bamtools_internal})
target_link_libraries(bamtools z)
