cmake_minimum_required(VERSION 3.9)
project(pydelt)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "6.4.0")
    message(FATAL_ERROR "Insufficient gcc version")
  endif()
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/pydeltser/pybind11)

pybind11_add_module(pydelt ./pydeltser/deltser_bindings.cpp)
target_include_directories(pydelt PRIVATE .)

if(MSVC)
    target_compile_options(pydelt PUBLIC $<$<CONFIG:RELEASE>: /Wall /O2>)
    target_compile_options(pydelt PUBLIC $<$<CONFIG:DEBUG>:/O1 /DEBUG:FULL /Zi /Zo>)
else()
    target_compile_options(pydelt PUBLIC $<$<CONFIG:RELEASE>: -pthread -O3>)
    target_compile_options(pydelt PUBLIC $<$<CONFIG:DEBUG>: -O2 -ggdb -D_GLIBCXX_DEBUG>)
endif()
