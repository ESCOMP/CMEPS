  class(client_type),                   intent(in) :: this  !< Pointer to the initialized client
  character(len=*),                     intent(in) :: key   !< The key to use to place the tensor
  integer, dimension(:),                intent(in) :: dims  !< Length along each dimension of the tensor