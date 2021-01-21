cmake_minimum_required(VERSION 3.17)

project(bouquet VERSION 1.0.0)

# link in the amrex lib
find_package(AMReX REQUIRED HINTS ${CMAKE_SOURCE_DIR}/libs/amrex)

# load mpi settings if needed
if($ENV{BOUQUET_MPI}==YES)
    find_package(MPI REQUIRED)

    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
endif()

# here we set special compiler flags
#   the cpp standard for this project is ++14
add_library(project_options INTERFACE)
target_compile_features(project_options INTERFACE cxx_std_17)

# include the all needed header files
include_directories(
        ./include
        ./libs/amrex/include
        ./src)

option(ENABLE_COVERAGE "Enable coverage reporting for  gcc/clang" FALSE)

if(ENABLE_COVERAGE)
    add_compile_options(--coverage -O0)
endif()

option(ENABLE_TESTING "Enable Test Builds" ON)

# include the dir of 3rd party libs 
# link_directories(/usr/lib64)
# link_directories(${CMAKE_SOURCE_DIR}/libs/amrex/lib)

add_subdirectory(./src)

# include test dir
if(ENABLE_TESTING)
    enable_testing()
    add_subdirectory(./tests)
endif()
