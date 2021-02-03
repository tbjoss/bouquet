
add_library(bouquet bouquet.cpp)

# link in amrex
target_link_libraries(bouquet AMReX::amrex AMReX::Flags_CXX)

add_executable(binary main.cpp)
target_link_libraries(binary PUBLIC project_options --coverage bouquet)