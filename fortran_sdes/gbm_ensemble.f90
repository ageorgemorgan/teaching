program solve_gbm 
    ! Simulates an ensemble of the geometric Brownian motion (GBM) following D. Higham 2001
    use stdlib_stats_distribution_normal, only: norm => rvs_normal
 
    implicit none
    real :: dt, rdt, winc, lambda = 2.0, mu = 1.0, x0 = 1.0, T = 1.0
    integer :: io, j, k, N_steps, N_steps_raw = 500, N_paths = 6, R = 4
    real, allocatable :: x(:, :), dw(:)

    dt = T/N_steps_raw
    N_steps = N_steps_raw/R
    rdt = R * dt ! called Dt in Higham's code but Fortran code is cAsE-inSENsiTive!
    
    allocate(dw(N_steps_raw))
    allocate(x(N_paths, N_steps))

    ! Initialize
    do j = 1, N_paths
        x(j, 1) = x0
    end do 

    !$omp parallel private(j, k, dw) shared(x) num_threads(6)

    !$omp do
    do j = 1, N_paths
        dw = sqrt(dt) * norm(0.0, 1.0, N_steps_raw)
        do k = 2, N_steps
            winc = sum(dw( R * (k - 1) + 1 : R * k))
            x(j, k) = x(j, k - 1) * ( 1.0 + rdt * lambda + mu * winc) 
        end do
    end do 
    !$omp end do 
    !$omp end parallel

    ! Record the state history
    open(newunit=io, file="gbm_sim.txt", status="replace", action="write")

    do j = 1, N_paths
        do k = 1, N_steps
            write(io, *) x(j, k)
        end do 
    end do
end program