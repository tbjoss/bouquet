add_executable(binary main.cpp)
target_link_libraries(binary PUBLIC project_options --coverage AMReX::amrex AMReX::Flags_CXX)