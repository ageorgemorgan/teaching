module time_stepper
   use, intrinsic :: iso_fortran_env, only: dp=>real64
   use read_utils
   implicit none
   private

   public :: do_time_stepping

contains
   subroutine do_time_stepping(initial_state_filename, history_filename, n_x, n_t, courant, method)
      ! Propagates an initial state on a spatial grid with n_x points out to n_t time-steps
      character(len=32), intent(inout) :: initial_state_filename ! inout bcz must go into a util subroutine
      integer, intent(inout) :: n_x, n_t ! inout bcz must go into a util subroutine
      character(len=32), intent(in) :: history_filename
      character(len=6), intent(in) :: method
      real(dp), intent(in) :: courant
      integer:: io, k_x, k_t
      real(dp), allocatable :: u0(:), u1(:)

      ! Load the initial state into the array u0
      call read_initial_state(initial_state_filename, u0, n_x)

      ! Allocate memory for u1
      allocate(u1(n_x))

      ! Create the outfile
      open(newunit=io, file=history_filename, status="replace", action="write")

      ! Record initial state
      do k_x = 1, n_x
         write(io, *) u0(k_x)
      end do

      ! March forward in time
      do k_t = 1, n_t - 1
         ! Get the new state

         ! CASE 1: Lax-Friedrichs Method
         if (method == "lf") then
            ! Left and right boundary nodes must be handled separately
            u1(1) = - 0.5_dp * courant * (u0(2)- u0(n_x)) + 0.5_dp * (u0(2) + u0(n_x))
            u1(n_x) = - 0.5_dp * courant * (u0(1) - u0(n_x - 1)) + 0.5_dp * (u0(1) + u0(n_x - 1))
            u1(2:n_x - 1) = - 0.5_dp * courant * (u0(3:n_x) - u0(1:n_x - 2)) + 0.5_dp * (u0(3:n_x) + u0(1:n_x-2))
         end if

         ! CASE 2: Lax-Wendroff Method
         if (method == "lw") then
            ! Left and right boundary nodes must be handled separately
            u1(1) = u0(1) - 0.5_dp * courant * (u0(2) - u0(n_x)) &
               + 0.5_dp * (courant**2) * (u0(2) - 2.0_dp * u0(1) + u0(n_x))
            u1(n_x) = u0(n_x) - 0.5_dp * courant * (u0(1) - u0(n_x - 1)) &
               + 0.5_dp * (courant**2) * (u0(1) - 2.0_dp * u0(n_x) + u0(n_x - 1))
            u1(2:n_x - 1) = u0(2:n_x - 1) - 0.5_dp * courant * (u0(3:n_x) - u0(1:n_x - 2)) &
               + 0.5_dp * (courant**2) * (u0(3:n_x) - 2.0_dp * u0(2:n_x - 1) + u0(1:n_x - 2))
         end if

         ! CASE 3: Upwind Method
         if (method == "upwind") then
            ! The left boundary node must be handled separately
            u1(1) = u0(1) - courant * (u0(1) - u0(n_x))
            u1(2:) = u0(2:) - courant * (u0(2:) - u0(1:n_x-1))
         end if

         ! CASE 4: Beam-Warming Method
         if (method == "bw") then
            ! The two leftmost nodes must be handled separately
            u1(1) = u0(1) - 0.5_dp * courant * (3.0_dp * u0(1) - 4.0_dp * u0(n_x) + u0(n_x - 1)) &
               + (0.5_dp * courant**2) * (u0(1) - 2.0_dp * u0(n_x) + u0(n_x - 1))
            u1(2) = u0(2) - 0.5_dp * courant * (3.0_dp * u0(2) - 4.0_dp * u0(1) + u0(n_x)) &
               + (0.5_dp * courant**2) * (u0(2) - 2.0_dp * u0(1) + u0(n_x))
            u1(3:) = u0(3:) - 0.5_dp * courant * (3.0_dp * u0(3:) - 4.0_dp * u0(2:n_x - 1) + u0(1:n_x - 2)) &
               + (0.5_dp * courant**2) * (u0(3:) - 2.0_dp * u0(2:n_x - 1) + u0(1:n_x - 2))
         end if

         ! Record the new state
         do k_x = 1, n_x
            write(io, *) u1(k_x)
         end do

         ! Set the old state equal to the new state
         u0 = u1
      end do

      ! Close the outfile
      close(io)
   end subroutine

end module
