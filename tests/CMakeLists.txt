add_subdirectory(./heat_equation_ex1)
add_subdirectory(./adv_amr_core)


#--------------------------------------------------------------------------------------------------
# hello world
#--------------------------------------------------------------------------------------------------
add_executable(000-test_hello_world 000-test-hello_world.cpp)
target_link_libraries(000-test_hello_world PUBLIC project_options --coverage bouquet)

add_test(000-tester 000-test_hello_world)

#--------------------------------------------------------------------------------------------------
# example grids
#--------------------------------------------------------------------------------------------------
add_executable(001-test_example_grid 001-test-example_grids.cpp)
target_link_libraries(001-test_example_grid PUBLIC project_options --coverage bouquet)

add_test(001-tester 001-test_example_grid)


