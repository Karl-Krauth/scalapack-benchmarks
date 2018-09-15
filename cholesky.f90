#define GAUSSIAN_A 1

program cholesky
    implicit none
    external blacs_exit
    external blacs_gridexit
    external blacs_gridinfo
    external descinit
    external sl_init
    external pdpotrf 
    external dlarnv
    external pdgeadd 
    integer, external :: numroc
    double precision, external :: MPI_Wtime 
    integer, parameter :: descriptor_len = 9
    double precision, parameter :: one = 1.0
    integer :: M 
    integer :: block_size 
    integer :: N_B 
    double precision :: lambda 
    integer :: processor_rows
    integer :: processor_cols 
    integer :: block_size_N_B
    integer :: context
    integer :: my_row
    integer :: my_col
    integer :: local_M
    integer :: local_N_A
    integer :: local_N_B
    integer :: leading_dim
    integer :: info
    integer :: descriptor_A(descriptor_len)
    integer :: descriptor_A_copy(descriptor_len)
    integer :: descriptor_B(descriptor_len)
    integer :: seed(4) = [0, 0, 0, 0]
    integer :: i
    double precision :: start_time
    double precision :: end_time
    double precision, allocatable :: temp_arr(:)
    double precision, allocatable :: A(:, :)
    double precision, allocatable :: A_copy(:, :)
    double precision, allocatable :: B(:, :)

    open(unit=1, file="out.txt")
    open(unit=2, file="in.txt")
    read (2, *) M, block_size, N_B, lambda, processor_rows, processor_cols
    ! Initialize the process grid.
    call sl_init(processor_rows, processor_cols, context)
    call blacs_gridinfo(context, processor_rows, processor_cols, my_row, my_col)

    ! Initialize the random number seed.
    seed(1) = mod(my_row, 4096)
    seed(2) = mod(my_col, 4096)
    seed(4) = 1

    ! Check if this process is on the process grid.
    if (my_row == -1) then
        go to 10
    end if

    if (my_row == 0 .and. my_col == 0) then
        write (1, *), "Num rows", M
        write (1, *), "Block size", block_size
        write (1, *), "B cols", N_B
        write (1, *), "lambda", lambda
        write (1, *), "processor rows", processor_rows
        write (1, *), "processor cols", processor_cols
    end if

    ! Compute matrix shapes.
    block_size_N_B = min(block_size, N_B)
    local_M = numroc(M, block_size, my_row, 0, processor_rows)
    local_N_A = numroc(M, block_size, my_col, 0, processor_cols)
    local_N_B = numroc(N_B, block_size, my_col, 0, processor_cols)
    leading_dim = max(1, local_M)

    ! Allocate local matrices.
    allocate(A(1:local_M, 1:local_N_A))
    allocate(B(1:local_M, 1:local_N_B))
#if GAUSSIAN_A
    allocate(A_copy(1:local_M, 1:local_N_A))
#endif

    ! Initialize global matrix descriptors.
    call descinit(descriptor_A, M, M, block_size, block_size, 0, 0, context, leading_dim, info)
    if (info /= 0) then
        write (1, *) "Descinit A failed argument", info, "is illegal."
        go to 10
    end if

    call descinit(descriptor_A_copy, M, M, block_size, block_size, 0, 0, context, leading_dim, info)
    if (info /= 0) then
        write(1, *) "Descinit A_copy failed argument", info, "is illegal."
        go to 10
    end if

    call descinit(descriptor_B, M, N_B, block_size, block_size_N_B, 0, 0, context, leading_dim, info)
    if (info /= 0) then
        write(1, *) "Descinit B failed argument", info, "is illegal."
        go to 10
    end if

    ! Set A = lambda I.
    if (my_row == 0 .and. my_col == 0) then
        print *, "Initializing A."
    endif
    call pdlaset("", M, M, 0.0, lambda, A, 1, 1, descriptor_A)

#if GAUSSIAN_A
    if (my_row == 0 .and. my_col == 0) then
        print *, "Adding noise to A."
    endif
    ! Add gaussian noise to all entries of A.
    allocate(temp_arr(1:local_N_A))
    do i = 1, local_M
        call dlarnv(3, seed, local_N_A, temp_arr) 
        A(i, :) = A(i, :) + temp_arr
    end do
    A_copy = A

    ! Set A = (A + A^t).
    call pdgeadd("T", M, M, one, A_copy, 1, 1, descriptor_A_copy, one, A, 1, 1, descriptor_A)
#endif


    ! Perform the cholesky.
    if (my_row == 0 .and. my_col == 0) then
        print *, "Running cholesky."
    endif
    start_time = MPI_Wtime()
    call pdpotrf("L", M, A, 1, 1, descriptor_A, info)
    end_time = MPI_Wtime()
    if (info /= 0) then
        write(1, *) "Cholesky failed leading minor of order", info, "is not positive definite."
        go to 10
    end if

    if (my_row == 0 .and. my_col == 0) then
        write(1, *) "Cholesky took", end_time - start_time, "seconds."
    end if

    ! Set B to be gaussian with mean 10.
    do i = 1, local_M
        call dlarnv(3, seed, local_N_B, B(i, :)) 
        B(i, :) = B(i, :) + 10.0
    end do

    10 continue
    call blacs_exit(0)
end program cholesky

subroutine sl_init(processor_rows, processor_cols, context)
    external blacs_get
    external blacs_gridinit
    external blacs_setup
    integer, intent(in) :: processor_cols
    integer, intent(in) :: processor_rows
    integer, intent(out) :: context
    integer my_id, num_procs

    call blacs_pinfo(my_id, num_procs)
    call blacs_get(-1, 0, context)
    call blacs_gridinit(context, 'Row-major', processor_rows, processor_cols)
    return
end subroutine sl_init
