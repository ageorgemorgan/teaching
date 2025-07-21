module read_utils
   use, intrinsic :: iso_fortran_env, only: dp=>real64, stderr => error_unit
   implicit none
   private

   ! Set up spacetime grid type
   type :: stgrid_type
      integer :: n_x, n_t
   end type

   ! Set up error info type
   type :: error_info_type
      integer :: nml_stat = 0, validation_stat = 0 ! status parameters: initially assume no errors occur
      character(len=512) :: nml_msg
      character(len=100), dimension(4) :: validation_msg = ["", "", "", ""] ! insert blank dummy messages
   contains
      procedure :: handle_nml_error, handle_validation_error
   end type

   public :: stgrid_type, error_info_type, read_namelist, read_initial_state

contains

   integer function handle_nml_error(self)
      ! Procedure for handling errors in finding the namelist
      class(error_info_type), intent(in) :: self ! note how we use class(...), and not type(...)

      if (self%nml_stat /= 0) then
         print *, "ERROR: "//self%nml_msg
         stop
      end if

      return

   end function

   integer function handle_validation_error(self)
      ! Procedure for handling validation errors in the namelist content
      ! These errors are: bad signal speed, bad Courant number, bad method,
      ! and bad grid parameters
      class(error_info_type), intent(in) :: self
      character(len=100), dimension(:), allocatable :: msg
      integer :: k

      if (self%validation_stat /= 0) then
         ! First throw out any "dummy" error messages
         msg = pack(self%validation_msg, self%validation_msg/="")

         do k = 1, size(msg)
            print *, "ERROR: "//msg(k)
         end do

         stop
      end if

      return

   end function

   subroutine read_namelist(file_path, stgrid, c, courant, initial_state_filename, history_filename, method, error_info)
      ! Reads Namelist from given file.
      ! Partially inspired by https://cyber.dabamos.de/programming/modernfortran/namelists.html
      character(len=*),  intent(in)    :: file_path
      type(stgrid_type), intent(inout) :: stgrid
      real(dp), intent(inout) :: c, courant
      character(len=32), intent(inout) :: initial_state_filename, history_filename
      character(len=6), intent(inout) :: method
      type(error_info_type), intent(out) :: error_info
      integer :: io, nml_stat
      character(len=512) :: nml_msg

      ! Namelist definition.
      namelist /config_sim/ stgrid, c, courant, initial_state_filename, history_filename, method

      ! Open and read Namelist file.
      open (action='read', file=file_path, iostat=nml_stat, newunit=io, iomsg=nml_msg)
      read (nml=config_sim, iostat=nml_stat, unit=io)

      if (nml_stat /= 0) then
         error_info%nml_stat = 1
         error_info%nml_msg = nml_msg
         close (io)
      end if

      ! validate Courant number
      if (courant < 0.0) then
         error_info%validation_stat = 1
         error_info%validation_msg(1) = "A negative Courant number is not allowed."
      end if
      if (courant > 1.0 .and. method/="bw") then
         error_info%validation_stat = 1
         error_info%validation_msg(1) = "Courant number is greater than 1, will cause numerical instability."
      end if
      ! Beam-Warming method is special since its maximum Courant number is 2
      if (courant > 2.0 .and. method=="bw") then
         error_info%validation_stat = 1
         error_info%validation_msg(1) = "Courant number is greater than 2, will cause numerical instability."
      end if

      ! Validate speed c
      if (c < 0.0) then
         error_info%validation_stat = 1
         error_info%validation_msg(2) = "A negative signal speed c is not allowed."
      end if

      ! Validate method keyword
      if (method /= "upwind" .and. method /= "lf" .and. method /= "lw" .and. method /= "bw") then
         error_info%validation_stat = 1
         error_info%validation_msg(3) = "Unrecognized method keyword. Acceptable values: 'upwind', 'lf', 'lw', 'bw'."
         ! Since "method" is only stored as a string of length 6, strictly speaking we aren't capturing the event where
         ! the user inputs "upwindasdf" or something like that: any chars appearing after the first six are automatically
         ! truncated. In my opinion this is not a big deal, and I'd rather be precise about storage of "method"
      end if

      ! Validate grid sizes
      if (stgrid%n_t <= 0 .or. stgrid%n_x <= 0) then
         error_info%validation_stat = 1
         error_info%validation_msg(4) = "Values of n_t, n_x must be strictly positive."
      end if


      close (io)
   end subroutine

   subroutine read_initial_state(initial_state_filename, u0, n_x)
      ! Reads an initial state and stores it in the array u0
      character(len=32), intent(in) :: initial_state_filename
      integer, intent(in) :: n_x
      real(dp), allocatable, intent(out) :: u0(:)
      integer :: io, k

      allocate(u0(n_x))

      open(newunit=io, file=initial_state_filename, action="read")

      do k = 1, n_x
         read(io, *) u0(k)
      end do

      close(io)

   end subroutine

end module
