#  Copyright (c) 2019, 2020, 2021 Tobias Brink
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#
# Required preamble
#

cmake_minimum_required(VERSION 3.10)

project(mytest)

find_package(KIM-API 2.2 REQUIRED CONFIG)

# End preamble

set(TEST_NAME "mytest")

add_executable(${TEST_NAME} core.cpp elastic.cpp dsl.cpp mytest.cpp
                            numer_forces.cpp poisson_disk_sampling.cpp)
set_target_properties(${TEST_NAME} PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)
target_link_libraries(${TEST_NAME} KIM-API::kim-api nlopt)


add_executable(test_switch_params core.cpp elastic.cpp dsl.cpp
                                  test_switch_params.cpp
                                  numer_forces.cpp poisson_disk_sampling.cpp)
set_target_properties(test_switch_params PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)
target_link_libraries(test_switch_params KIM-API::kim-api nlopt)



add_executable(test_unit_conv core.cpp elastic.cpp dsl.cpp
                              test_unit_conv.cpp
                              numer_forces.cpp poisson_disk_sampling.cpp)
set_target_properties(test_unit_conv PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)
target_link_libraries(test_unit_conv KIM-API::kim-api nlopt)


add_executable(test_benchmark core.cpp elastic.cpp dsl.cpp
                              test_benchmark.cpp
                              numer_forces.cpp poisson_disk_sampling.cpp)
set_target_properties(test_benchmark PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)
target_link_libraries(test_benchmark KIM-API::kim-api nlopt)


add_executable(test_write_params core.cpp elastic.cpp dsl.cpp
                                 test_write_params.cpp
                                 numer_forces.cpp poisson_disk_sampling.cpp)
set_target_properties(test_write_params PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED ON
)
target_link_libraries(test_write_params KIM-API::kim-api nlopt)
