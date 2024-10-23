cmake_minimum_required(VERSION 3.15)

project(ambi~)

set(PDCMAKE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/Resources/pd.cmake/"
    CACHE PATH "Path to pd.cmake")
include(${PDCMAKE_DIR}/pd.cmake)

# ╭──────────────────────────────────────╮
# │               LibSofa                │
# ╰──────────────────────────────────────╯

cmake_policy(SET CMP0077 NEW)
set(BUILD_TESTS OFF)
option(BUILD_SHARED_LIBS "Build shared library" OFF)

set(BUILD_STATIC_LIBS ON)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libmysofa/
                 EXCLUDE_FROM_ALL)
set_target_properties(mysofa-static PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Set variables for libmysofa
set(MYSOFA_INCLUDE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libmysofa/src/hrtf")
set(MYSOFA_LIBRARIES mysofa-static)
include_directories(${MYSOFA_INCLUDE_DIRS})

# ╭──────────────────────────────────────╮
# │           spatialaudioLib            │
# ╰──────────────────────────────────────╯
set(BUILD_SHARED_LIBS OFF)
set(MYSOFA_FOUND ON)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio
                 EXCLUDE_FROM_ALL)
set_target_properties(spatialaudio-static PROPERTIES POSITION_INDEPENDENT_CODE
                                                     ON)

# ╭──────────────────────────────────────╮
# │             Pd External              │
# ╰──────────────────────────────────────╯
file(GLOB EXTERNAL_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Sources/*.cpp")
pd_add_external(ambi~ ${EXTERNAL_SOURCES} TARGET ambi_tilde)

add_definitions(-DHAVE_MYSOFA)
target_include_directories(
  ambi_tilde
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/include/hrtf/")
target_include_directories(
  ambi_tilde
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/include/")
target_include_directories(
  ambi_tilde
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/include/dsp")
target_include_directories(
  ambi_tilde
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/source/kiss_fft/")

target_include_directories(
  ambi_tilde PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/Libraries/libmysofa/src/")
set_target_properties(ambi_tilde PROPERTIES CXX_STANDARD 17)
target_link_libraries(ambi_tilde PRIVATE spatialaudio-static)

# ╭──────────────────────────────────────╮
# │              Data Files              │
# ╰──────────────────────────────────────╯
pd_add_datafile(ambi_tilde
                "${CMAKE_CURRENT_SOURCE_DIR}/Resources/ambi~-help.pd")
pd_add_datafile(ambi_tilde "${CMAKE_CURRENT_SOURCE_DIR}/Resources/nh906.sofa")
