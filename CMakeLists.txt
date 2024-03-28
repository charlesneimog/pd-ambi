# Define your standard CMake header (for example):
cmake_minimum_required(VERSION 3.5)
include(${CMAKE_CURRENT_SOURCE_DIR}/Resources/pd.build/pd.cmake)

project(ambi)

set_pd_external_path(${CMAKE_CURRENT_SOURCE_DIR}/ambi~)

set(PD_FLOATSIZE64 "OFF")

#╭──────────────────────────────────────╮
#│               Kiss FFT               │
#╰──────────────────────────────────────╯
file(GLOB KISSFFT_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/source/kiss_fft/*.c")
add_library(kissfft STATIC ${KISSFFT_SOURCES})

#╭──────────────────────────────────────╮
#│           spatialaudioLib            │
#╰──────────────────────────────────────╯
cmake_policy(SET CMP0077 NEW)
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio)

#╭──────────────────────────────────────╮
#│             Pd External              │
#╰──────────────────────────────────────╯
file(GLOB EXTERNAL_SOURCES "${CMAKE_SOURCE_DIR}/Sources/*.cpp")
add_pd_external(elseambi ambi~ ${EXTERNAL_SOURCES})
target_include_directories(elseambi PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/include")
target_include_directories(elseambi PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/Libraries/libspatialaudio/source/kiss_fft/")
target_link_libraries(elseambi spatialaudio-static kissfft)

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
