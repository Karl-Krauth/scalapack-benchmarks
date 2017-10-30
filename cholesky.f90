program cholesky
external blacs_exit
external blacs_gridexit
external blacs_gridinfo
external descinit
external matinit
external sl_init
integer, parameter :: N = 16384
integer, parameter :: block_size = 4096
integer :: processor_rows = 4
integer :: processor_cols = 1
integer :: context
integer :: my_row
integer :: my_col

! Initialize the process grid.
call sl_init(processor_rows, processor_cols, context)
call blacs_gridinfo(context, processor_rows, processor_cols, my_row, my_col)

! Check if this process is on the process grid.
if (my_row .eq. -1) then
    go to 10
end if

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
