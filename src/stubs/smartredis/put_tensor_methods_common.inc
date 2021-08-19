  class(client_type),                    intent(in) :: this !< Fortran SLIC client
  character(len=*),                      intent(in) :: key  !< The unique key used to store in the database
  integer, dimension(:),                 intent(in) :: dims !< The length of each dimension
  