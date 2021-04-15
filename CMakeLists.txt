cmake_minimum_required(VERSION 3.0)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

project (pfbwt-f)

find_library(HTS_LIB hts)

# Google Test support
set(CMAKE_CXX_STANDARD 17)
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/6c5c4554ac218a8e19168edc121b1ad232015185.zip
)
FetchContent_MakeAvailable(googletest)
include(CTest)
enable_testing()

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/sdsl-lite/include)

message("Building in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_FLAGS "-DM64  -Wall -Wextra -std=c++17")
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

# tests
add_executable(marker_test tests/marker_test.cpp)
target_link_libraries(marker_test gtest_main)
add_test(NAME single_chrom_no_markers COMMAND bash ${PROJECT_SOURCE_DIR}/tests/vcf_to_bwt_test_no_markers.sh ${PROJECT_SOURCE_DIR} single_chrom)
add_test(NAME mult_chroms_no_markers COMMAND bash ${PROJECT_SOURCE_DIR}/tests/vcf_to_bwt_test_no_markers.sh ${PROJECT_SOURCE_DIR} mult_chroms)
add_test(NAME single_chrom COMMAND bash ${PROJECT_SOURCE_DIR}/tests/vcf_to_bwt_test.sh ${PROJECT_SOURCE_DIR} single_chrom)
add_test(NAME mult_chroms COMMAND bash ${PROJECT_SOURCE_DIR}/tests/vcf_to_bwt_test.sh ${PROJECT_SOURCE_DIR} mult_chroms)
include(GoogleTest)
gtest_discover_tests(marker_test)
include(CPack)

# file(COPY ${PROJECT_SOURCE_DIR}/vcf_to_bwt.py DESTINATION bin FILE_PERMISSIONS OWNER_READ OWNER_EXECUTE)
install(FILES ${PROJECT_SOURCE_DIR}/vcf_to_bwt.py DESTINATION bin PERMISSIONS OWNER_READ OWNER_EXECUTE)
install(TARGETS pfbwt-f64 merge_pfp vcf_scan merge_mps dump_markers mps_to_ma
        DESTINATION bin)
