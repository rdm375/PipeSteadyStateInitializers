cmake_minimum_required(VERSION 3.16)
project(isothermal_steady_initializers LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(pipe_initializers
    src/steady_initializers.cpp
    src/adiabatic_initializers.cpp
    src/thermal_wall_initializers.cpp
    src/soil_shell_initializers.cpp)
target_include_directories(pipe_initializers PUBLIC include)
target_compile_options(pipe_initializers PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(pipe_init_profile_export examples/export_initializer_profiles.cpp)
target_link_libraries(pipe_init_profile_export PRIVATE pipe_initializers)
target_compile_options(pipe_init_profile_export PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(pipe_init_demo src/main.cpp)
target_link_libraries(pipe_init_demo PRIVATE pipe_initializers)
target_compile_options(pipe_init_demo PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(adiabatic_pipe_init_demo src/adiabatic_demo.cpp)
target_link_libraries(adiabatic_pipe_init_demo PRIVATE pipe_initializers)
target_compile_options(adiabatic_pipe_init_demo PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(adiabatic_pipe_init_profile_export examples/export_adiabatic_profiles.cpp)
target_link_libraries(adiabatic_pipe_init_profile_export PRIVATE pipe_initializers)
target_compile_options(adiabatic_pipe_init_profile_export PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(thermal_wall_init_demo src/thermal_wall_demo.cpp)
target_link_libraries(thermal_wall_init_demo PRIVATE pipe_initializers)
target_compile_options(thermal_wall_init_demo PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(thermal_wall_profile_export examples/export_thermal_wall_profiles.cpp)
target_link_libraries(thermal_wall_profile_export PRIVATE pipe_initializers)
target_compile_options(thermal_wall_profile_export PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(soil_shell_init_demo src/soil_shell_demo.cpp)
target_link_libraries(soil_shell_init_demo PRIVATE pipe_initializers)
target_compile_options(soil_shell_init_demo PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(soil_shell_profile_export examples/export_soil_shell_profiles.cpp)
target_link_libraries(soil_shell_profile_export PRIVATE pipe_initializers)
target_compile_options(soil_shell_profile_export PRIVATE -Wall -Wextra -Wpedantic -Werror)

add_executable(pipe_init_benchmark examples/benchmark_initializers.cpp)
target_link_libraries(pipe_init_benchmark PRIVATE pipe_initializers)
target_compile_options(pipe_init_benchmark PRIVATE -Wall -Wextra -Wpedantic -Werror)

enable_testing()
add_executable(pipe_init_tests tests/test_initializers.cpp)
target_link_libraries(pipe_init_tests PRIVATE pipe_initializers)
target_compile_options(pipe_init_tests PRIVATE -Wall -Wextra -Wpedantic -Werror)
add_test(NAME pipe_init_tests COMMAND pipe_init_tests)
