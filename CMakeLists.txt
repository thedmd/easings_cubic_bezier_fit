cmake_minimum_required(VERSION 3.16)

project(bezier_fit VERSION 0.1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(SOURCES
    src/main.cpp
    src/FitCurves.hpp
    src/FitCurves.c
    src/bezier_fit.natvis
)

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})

add_executable(${PROJECT_NAME} ${SOURCES})

find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PRIVATE
        OpenMP::OpenMP_CXX
    )
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin
)

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT bezier_fit)

