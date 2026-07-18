cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(twoCellShear LANGUAGES C CXX)

if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(
    HEMOCELL_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../src_hemocell"
    CACHE PATH
    "Path to the HemoCell source directory"
)
get_filename_component(HEMOCELL_SOURCE_DIR "${HEMOCELL_SOURCE_DIR}" ABSOLUTE)

set(HEMOCELL_LIBRARY "${HEMOCELL_SOURCE_DIR}/build/libhemocell.a")
if(NOT EXISTS "${HEMOCELL_LIBRARY}")
    message(FATAL_ERROR
        "HemoCell library not found at ${HEMOCELL_LIBRARY}. "
        "Run src_hemocell/compile.sh first."
    )
endif()

set(HDF5_PREFER_PARALLEL TRUE)
find_package(MPI REQUIRED COMPONENTS CXX)
find_package(HDF5 REQUIRED COMPONENTS C HL)

add_library(hemocell::hemocell STATIC IMPORTED GLOBAL)
set_target_properties(hemocell::hemocell PROPERTIES
    IMPORTED_LOCATION "${HEMOCELL_LIBRARY}"
)

add_executable(twoCellShear twoCellShear.cpp)

target_compile_definitions(twoCellShear PRIVATE
    PLB_MPI_PARALLEL
    PLB_USE_POSIX
    PLB_SMP_PARALLEL
    OMPI_SKIP_MPICXX
)

target_include_directories(twoCellShear PRIVATE
    "${HEMOCELL_SOURCE_DIR}"
    "${HEMOCELL_SOURCE_DIR}/config"
    "${HEMOCELL_SOURCE_DIR}/core"
    "${HEMOCELL_SOURCE_DIR}/helper"
    "${HEMOCELL_SOURCE_DIR}/io"
    "${HEMOCELL_SOURCE_DIR}/mechanics"
    "${HEMOCELL_SOURCE_DIR}/external/tinyxml2"
    "${HEMOCELL_SOURCE_DIR}/palabos/src"
    "${HEMOCELL_SOURCE_DIR}/palabos/externalLibraries"
    "${HEMOCELL_SOURCE_DIR}/palabos/externalLibraries/tinyxml"
    ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(twoCellShear PRIVATE
    hemocell::hemocell
    MPI::MPI_CXX
    ${HDF5_HL_LIBRARIES}
    ${HDF5_LIBRARIES}
)

set_target_properties(twoCellShear PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
