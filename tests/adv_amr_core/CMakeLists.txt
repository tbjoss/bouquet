include_directories(./Src_K)
#--------------------------------------------------------------------------------------------------
# example advection amr core
#--------------------------------------------------------------------------------------------------
add_executable(example_adv_amr_core AdvancePhiAllLevels.cpp AdvancePhiAtLevel.cpp AmrCoreAdv.cpp DefineVelocity.cpp main.cpp)
target_link_libraries(example_adv_amr_core PUBLIC project_options --coverage bouquet)

add_test(tester_adv_amr_core example_adv_amr_core "${CMAKE_SOURCE_DIR}/tests/adv_amr_core/inputs")