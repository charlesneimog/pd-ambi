cmake_minimum_required(VERSION 3.15)
include(${CMAKE_CURRENT_SOURCE_DIR}/Resources/pd.build/pd.cmake)

project(ambi)
set_pd_external_path(${CMAKE_CURRENT_SOURCE_DIR}/ambi~)
set(PD_FLOATSIZE64 "OFF")

#╭──────────────────────────────────────╮
#│               LibSofa                │
#╰──────────────────────────────────────╯
cmake_policy(SET CMP0077 NEW)
set(BUILD_TESTS OFF)
option(BUILD_SHARED_LIBS "Build shared library" OFF)

set(BUILD_STATIC_LIBS ON)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libmysofa/)
set_target_properties(mysofa-static PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Set variables for libmysofa
set(MYSOFA_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libmysofa/src/hrtf")
set(MYSOFA_LIBRARIES mysofa-static)
include_directories(${MYSOFA_INCLUDE_DIRS})

#╭──────────────────────────────────────╮
#│           spatialaudioLib            │
#╰──────────────────────────────────────╯
set(BUILD_SHARED_LIBS OFF)
set(MYSOFA_FOUND ON)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio)
set_target_properties(spatialaudio-static PROPERTIES POSITION_INDEPENDENT_CODE ON)

#╭──────────────────────────────────────╮
#│             Pd External              │
#╰──────────────────────────────────────╯
file(GLOB EXTERNAL_SOURCES "${CMAKE_SOURCE_DIR}/Sources/*.cpp")
add_pd_external(ambi_tilde ambi~ ${EXTERNAL_SOURCES})

add_definitions(-DHAVE_MYSOFA)
target_include_directories(ambi_tilde PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/include/hrtf/")
target_include_directories(ambi_tilde PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/include/")
target_include_directories(ambi_tilde PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/source/kiss_fft/")
target_include_directories(ambi_tilde PRIVATE "${CMAKE_BINARY_DIR}/Libraries/libmysofa/src/")
set_target_properties(ambi_tilde PROPERTIES CXX_STANDARD 11)
target_link_libraries(ambi_tilde PRIVATE spatialaudio-static)

#╭──────────────────────────────────────╮
#│    Copy Help Patch and .sofa file    │
#╰──────────────────────────────────────╯


if (LINUX)
    add_custom_target(Extra_Files ALL
            COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/Resources/ambi~-help.pd" "${CMAKE_CURRENT_SOURCE_DIR}/ambi~/"
            COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/Resources/nh906.sofa" "${CMAKE_CURRENT_SOURCE_DIR}/ambi~/"
            DEPENDS ambi_tilde)
endif()





