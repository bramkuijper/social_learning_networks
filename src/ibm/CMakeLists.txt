# this file is used for the programme CMake,
# which allows us to build and compile the executable
# see https://cmake.org/cmake/help/latest/index.html
cmake_minimum_required(VERSION 3.10)

set(CMAKE_BUILD_TYPE Debug)

#### VARIABLES ####

# put executable name in variable, as we
# need it in multiple places
set(EXE_NAME "NetworksIBM")


# google test url
set(GOOGLETEST_URL "https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip") 

#### THE ACTUAL PROJECT FILE ####

# set the project name and version
project(NetworksIBM)

add_executable("${EXE_NAME}" network_main.cpp)
add_library(network_library network_main.cpp patch.cpp individual.cpp simulation.cpp)

# add the current dir in which the executable is being built
# as include dir, as there are a bunch of header files in there
target_include_directories(
    "${EXE_NAME}" 
    PUBLIC
    "${PROJECT_BINARY_DIR}"
    )

# make sure it is c++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

# now adding the google mock testing framework
# see: https://google.github.io/googletest/quickstart-cmake.html

# call cmake's FetchContent module
# see https://cmake.org/cmake/help/latest/module/FetchContent.html 
include(FetchContent)
FetchContent_Declare(
    googletest
    URL "${GOOGLETEST_URL}"
)

FetchContent_MakeAvailable(googletest)

enable_testing()

# link libraries to main exe
target_link_libraries(
    "${EXE_NAME}"
    network_library
)

target_include_directories(
    "${EXE_NAME}" 
    PUBLIC
    "${PROJECT_BINARY_DIR}"
    )

include(GoogleTest)

