cmake_minimum_required(VERSION 2.8.12)
project(pydelt)
add_subdirectory(./pydeltser/pybind11)
pybind11_add_module(pydelt ./pydeltser/deltser_bindings.cpp)
