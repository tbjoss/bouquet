#--------------------------------------------------------------------------------------------------
# example heat equation
#--------------------------------------------------------------------------------------------------
add_executable(heat_equation myfunc.cpp main.cpp)
target_link_libraries(heat_equation PUBLIC project_options --coverage bouquet)

add_test(tester-heat_equation heat_equation "${CMAKE_SOURCE_DIR}/tests/heat_equation_ex1/inputs")