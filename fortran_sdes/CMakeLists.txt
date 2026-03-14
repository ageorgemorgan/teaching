cmake_minimum_required(VERSION 3.10)

find_package(fortran_stdlib REQUIRED)

set(PROJECT_NAME baby_sde)

project(baby_sde Fortran)
set(WEINER_PROGRAM_NAME wiener_ensemble.exe)
set(GBM_PROGRAM_NAME gbm_ensemble.exe)
add_executable(${WEINER_PROGRAM_NAME} wiener_ensemble.f90)
add_executable(${GBM_PROGRAM_NAME} gbm_ensemble.f90)

target_link_libraries(
  ${WEINER_PROGRAM_NAME}
  PRIVATE
  fortran_stdlib::fortran_stdlib
)

target_link_libraries(
  ${GBM_PROGRAM_NAME}
  PRIVATE
  fortran_stdlib::fortran_stdlib
)
