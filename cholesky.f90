program cholesky
implicit none
external blacs_exit
external blacs_gridexit
external blacs_gridinfo
external descinit
external matinit
external sl_init
integer, external :: numroc
integer, parameter :: M = 16384
integer, parameter :: block_size = 4096
integer, parameter :: N_B = 1
integer, parameter :: descriptor_len = 9
integer :: block_size_N_B
integer :: processor_rows = 4
integer :: processor_cols = 1
integer :: context
integer :: my_row
integer :: my_col
integer :: local_M
integer :: local_N_A
integer :: local_N_B
integer :: leading_dim
integer :: info
integer :: descriptor_A(descriptor_len)
integer :: descriptor_B(descriptor_len)
integer :: seed(4) = [0, 0, 0, 0]
double precision, allocatable :: A(:, :)
double precision, allocatable :: B(:, :)

! Initialize the process grid.
call sl_init(processor_rows, processor_cols, context)
call blacs_gridinfo(context, processor_rows, processor_cols, my_row, my_col)

! Initialize the random number seed.
seed(1) = my_row
seed(2) = my_col

! Check if this process is on the process grid.
if (my_row == -1) then
    go to 10
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

! Initialize global matrix descriptors.
call descinit(descriptor_A, M, M, block_size, block_size, 0, 0, context, leading_dim, info)
if (info /= 0) then
    print *, "Descinit A failed argument", info, "is illegal."
    go to 10
end if

call descinit(descriptor_B, M, N_B, block_size, block_size_N_B, 0, 0, context, leading_dim, info)
if (info /= 0) then
    print *, "Descinit B failed argument", info, "is illegal."
    go to 10
end if

! Randomly initialize local matrices.
call dlarnv(3, seed, local_M * local_N_A, A)
call dlarnv(3, seed, local_M * local_N_B, B)

! Compute A = A + A^t
call pdgeadd("T", M, M, 1.0, A, 1, 1, descriptor_A, 1.0, A, 1, 1, descriptor_A)
10 continue
call blacs_exit(0)
write (*, *), "DONE!"

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
