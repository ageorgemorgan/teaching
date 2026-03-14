cmake_minimum_required(VERSION 3.10)

find_package(fortran_stdlib REQUIRED)

set(PROJECT_NAME baby_sde)

project(baby_sde Fortran)
set(PROGRAM_NAME baby_sde.exe)
add_executable(${PROGRAM_NAME} wiener_ensemble.f90)

target_link_libraries(
  ${PROGRAM_NAME}
  PRIVATE
  fortran_stdlib::fortran_stdlib
)
