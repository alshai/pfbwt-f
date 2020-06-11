cmake_minimum_required(VERSION 3.0)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

project (pfbwt-f)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/sdsl-lite/include)

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "-DM64  -Wall -Wextra -std=c++11")
set(CMAKE_C_FLAGS "-DM64 -Wall -Wextra -std=c99")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-g -Ofast -fstrict-aliasing -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")

add_executable(pfbwt-f64 src/pfbwt-f.cpp gsa/gsacak.c src/utils.c)
TARGET_LINK_LIBRARIES(pfbwt-f64 z)
add_executable(merge_pfp src/merge_pfp.cpp gsa/gsacak.c src/utils.c)
TARGET_LINK_LIBRARIES(merge_pfp z pthread)
add_executable(vcf_scan src/vcf_scan.cpp)
TARGET_LINK_LIBRARIES(vcf_scan hts)
add_executable(merge_mps src/merge_mps.cpp)
add_executable(dump_markers src/dump_markers.cpp)
add_executable(mps_to_ma src/mps_to_ma.cpp src/utils.c)

# file(COPY ${PROJECT_SOURCE_DIR}/vcf_to_bwt.py DESTINATION bin FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)
install(FILES ${PROJECT_SOURCE_DIR}/vcf_to_bwt.py DESTINATION bin PERMISSIONS OWNER_READ OWNER_EXECUTE)
install(TARGETS pfbwt-f64 merge_pfp vcf_scan merge_mps dump_markers mps_to_ma
        DESTINATION bin)
