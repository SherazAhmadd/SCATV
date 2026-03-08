cmake_minimum_required(VERSION 3.10)


project(SCATV_Aligner)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(SCATV 
    src/main.cpp 
    src/parser.cpp 
    src/index.cpp 
    src/align.cpp
)