# Define your standard CMake header (for example):
cmake_minimum_required(VERSION 3.15)
include(${CMAKE_CURRENT_SOURCE_DIR}/Resources/pd.build/pd.cmake)
project(ambi)

set_pd_external_path(${CMAKE_CURRENT_SOURCE_DIR}/ambi~)

set(PD_FLOATSIZE64 "OFF")

#╭──────────────────────────────────────╮
#│               LibSofa                │
#╰──────────────────────────────────────╯
cmake_policy(SET CMP0077 NEW)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libmysofa/)
set(BUILD_TESTS OFF)
set(BUILD_SHARED_LIBS OFF)
set(BUILD_STATIC_LIBS ON)
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
# include_directories(elseambi PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/include/hrtf/")
# include_directories(elseambi PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/source/kiss_fft/")


#╭──────────────────────────────────────╮
#│             Pd External              │
#╰──────────────────────────────────────╯
file(GLOB EXTERNAL_SOURCES "${CMAKE_SOURCE_DIR}/Sources/*.cpp")
add_pd_external(elseambi ambi~ ${EXTERNAL_SOURCES})
target_include_directories(elseambi PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/include/hrtf/")
target_include_directories(elseambi PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/include/")
target_include_directories(elseambi PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/source/kiss_fft/")
target_include_directories(elseambi PRIVATE "${CMAKE_BINARY_DIR}/Libraries/libmysofa/src/")
add_definitions(-DHAVE_MYSOFA)
target_link_libraries(elseambi spatialaudio-static)

# set_target_properties(elseambi PROPERTIES CXX_STANDARD 11)
# target_compile_options(elseambi PRIVATE 
#                 "-O3" 
#                 "-funroll-loops" 
#                 "-fomit-frame-pointer" 
#                 "-mfpmath=sse" 
#                 "-msse" 
#                 "-msse2" 
#                 "-msse3" 
#                 "-DPD"
# )
