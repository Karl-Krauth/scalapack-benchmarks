#define GAUSSIAN_A 1

program dgemm
    implicit none
    external blacs_exit
    external blacs_gridexit
    external blacs_gridinfo
    external blacs_barrier
    external descinit
    external sl_init
    external pdpotrf 
    external dlarnv
    external pdgeadd 
    integer, external :: indxg2p
    integer, external :: numroc
    double precision, external :: MPI_Wtime 
    integer, parameter :: descriptor_len = 9
    double precision, parameter :: one = 1.0
    integer :: M 
    integer :: block_size 
    integer :: processor_rows
    integer :: processor_cols 
    integer :: context
    integer :: my_row
    integer :: my_col
    integer :: local_M
    integer :: local_N
    integer :: work_size
    integer :: leading_dim
    integer :: info
    integer :: descriptor_A(descriptor_len)
    integer :: descriptor_B(descriptor_len)
    integer :: descriptor_C(descriptor_len)
    integer :: seed(4) = [0, 0, 0, 0]
    integer :: i
    double precision :: start_time
    double precision :: end_time
    double precision, allocatable :: A(:, :)
    double precision, allocatable :: B(:, :)
    double precision, allocatable :: C(:, :)

    open(unit=1, file="out.txt")
    open(unit=2, file="in.txt")
    read (2, *) M, block_size, processor_rows, processor_cols
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
        write (1, *), "processor rows", processor_rows
        write (1, *), "processor cols", processor_cols
    end if

    ! Compute matrix shapes.
    local_M = numroc(M, block_size, my_row, 0, processor_rows)
    local_N = numroc(M, block_size, my_col, 0, processor_cols)
    leading_dim = max(1, local_M)

    ! Allocate local matrices.
    allocate(A(1:local_M, 1:local_N))
    allocate(B(1:local_M, 1:local_N))
    allocate(C(1:local_M, 1:local_N))

    ! Initialize global matrix descriptors.
    call descinit(descriptor_A, M, M, block_size, block_size, 0, 0, context, leading_dim, info)
    if (info /= 0) then
        write (1, *) "Descinit A failed argument", info, "is illegal."
        go to 10
    end if

    call descinit(descriptor_B, M, M, block_size, block_size, 0, 0, context, leading_dim, info)
    if (info /= 0) then
        write (1, *) "Descinit B failed argument", info, "is illegal."
        go to 10
    end if

    ! Initialize global matrix descriptors.
    call descinit(descriptor_C, M, M, block_size, block_size, 0, 0, context, leading_dim, info)
    if (info /= 0) then
        write (1, *) "Descinit C failed argument", info, "is illegal."
        go to 10
    end if

    if (my_row == 0 .and. my_col == 0) then
        print *, "Adding noise to A and B."
    endif

    ! Add gaussian noise to all entries of A and B.
    call dlarnv(3, seed, local_M * local_N, A) 
    call dlarnv(3, seed, local_M * local_N, B) 

    ! Perform the DGEMM.
    if (my_row == 0 .and. my_col == 0) then
        print *, "Running DGEMM."
    endif
    start_time = MPI_Wtime()
    call pdgemm("N", "N", M, M, M, one, A, 1, 1, descriptor_A, B, 1, 1, descriptor_B, 0.0, C, 1, 1, descriptor_C)
    if (info /= 0) then
        write(1, *) "DGEMM failed with error code:", info
        go to 10
    end if
    call blacs_barrier(context, "A")
    end_time = MPI_Wtime()

    if (my_row == 0 .and. my_col == 0) then
        write(1, *) "DGEMM took", end_time - start_time, "seconds."
    end if

    10 continue
    call blacs_exit(0)
end program dgemm

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
