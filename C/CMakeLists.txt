cmake_minimum_required(VERSION 3.27)
project(HW1_C C)

set(CMAKE_C_STANDARD 11)

add_library(HW1_C STATIC
        kmeans.c)
