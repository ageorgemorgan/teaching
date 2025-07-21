program transport_solver
   use, intrinsic :: iso_fortran_env, only: dp=>real64
   use read_utils
   use time_stepper
   implicit none

   ! Simulation parameters to be read in from namelist
   type(stgrid_type) :: stgrid ! defined in read_utils
   real(dp) :: c, courant
   character(len=32) :: initial_state_filename, history_filename
   character(len=6) :: method
   type(error_info_type) :: error_info
   integer :: err_flag

   ! Load and validate the namelist
   call read_namelist("config_sim.nml", stgrid, c, courant, initial_state_filename, history_filename, method, error_info)

   ! Handle any errors w/ loading the namelist
   err_flag = error_info%handle_nml_error()

   ! Handle validation errors
   err_flag = error_info%handle_validation_error()

   ! Do time-stepping
   call do_time_stepping(initial_state_filename, history_filename, stgrid%n_x, stgrid%n_t, courant, method)
end program
