program wiener_sim
    ! Simulates an ensemble of the Wiener process following D. Higham 2001
    use stdlib_stats_distribution_normal, only: norm => rvs_normal
 
    implicit none
    integer :: io, j, k, N_steps = 500, N_paths = 6
    real :: dt, T = 1.0
    real, allocatable :: w(:, :), dw(:)

    dt = T/N_steps
    
    allocate(dw(N_steps))
    allocate(w(N_paths, N_steps))

    !$omp parallel private(k, dw) shared(w) num_threads(6)

    !$omp do
    do k = 1, N_paths
        dw = sqrt(dt) * norm(0.0, 1.0, N_steps)
        w(k, :) = cumsum(dw)
    end do 
    !$omp end do 
    !$omp end parallel

    ! Record the state history
    open(newunit=io, file="wiener_sim.txt", status="replace", action="write")

    do j = 1, N_paths
        do k = 1, N_steps
            write(io, *) w(j, k)
        end do 
    end do

contains

    function cumsum(x) result(csx)
        implicit none
        real, intent(in) :: x(:)
        integer :: n, k
        real, allocatable :: csx(:)

        n = size(x)

        allocate(csx(n))
      
        csx(1) = x(1)

        do k = 2, n
            csx(k) = csx(k - 1) + x(k)
        end do 

    end function cumsum

end program