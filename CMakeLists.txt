cmake_minimum_required(VERSION 3.10)
project(CentralityEngine)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(APPLE)
  if(IS_DIRECTORY /opt/homebrew/opt/llvm/bin)
    set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
    set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
  else()
    message(FATAL_ERROR "LLVM not found. Please run 'brew install llvm'")
  endif()
endif()

find_package(OpenMP REQUIRED)

add_executable(centrality_engine main.cpp centrality_engine.cpp)

target_link_libraries(centrality_engine PRIVATE OpenMP::OpenMP_CXX)

install(TARGETS centrality_engine DESTINATION bin)