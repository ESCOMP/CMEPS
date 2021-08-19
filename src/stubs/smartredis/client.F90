module smartredis_client

use iso_c_binding, only : c_ptr, c_bool, c_null_ptr, c_char, c_int
use iso_c_binding, only : c_int8_t, c_int16_t, c_int32_t, c_int64_t, c_float, c_double, c_size_t

implicit none; private

!> Stores all data and methods associated with the SmartRedis client that is used to communicate with the database
type, public :: client_type
  private

  contains

  ! Public procedures
  !> Puts a tensor into the database (overloaded)
  generic :: put_tensor => put_tensor_i8, put_tensor_i16, put_tensor_i32, put_tensor_i64, &
                           put_tensor_float, put_tensor_double
  !> Retrieve the tensor in the database into already allocated memory (overloaded)
  generic :: unpack_tensor => unpack_tensor_i8, unpack_tensor_i16, unpack_tensor_i32, unpack_tensor_i64, &
                              unpack_tensor_float, unpack_tensor_double

  !> Initializes a new instance of the SmartRedis client
  procedure :: initialize
  !> Initializes a new instance of the SmartRedis client
  procedure :: isinitialized
  
  !> Destructs a new instance of the SmartRedis client
  procedure :: destructor
  !> Check the database for the existence of a specific model
  procedure :: model_exists
  !> Check the database for the existence of a specific tensor
  procedure :: tensor_exists
  !> Check the database for the existence of a specific key
  procedure :: key_exists
  !> Poll the database and return if the model exists
  procedure :: poll_model
  !> Poll the database and return if the tensor exists
  procedure :: poll_tensor
  !> Poll the database and return if the key exists
  procedure :: poll_key
  !> Rename a tensor within the database
  procedure :: rename_tensor
  !> Delete a tensor from the database
  procedure :: delete_tensor
  !> Copy a tensor within the database to a new key
  procedure :: copy_tensor
  !> Set a model from a file
  procedure :: set_model_from_file
  !> Set a model from a byte string that has been loaded within the application
  procedure :: set_model
  !> Retrieve the model as a byte string
  procedure :: get_model
  !> Set a script from a specified file
  procedure :: set_script_from_file
  !> Set a script as a byte or text string
  procedure :: set_script
  !> Retrieve the script from the database
  procedure :: get_script
  !> Run a script that has already been stored in the database
  procedure :: run_script
  !> Run a model that has already been stored in the database
  procedure :: run_model
  !> Put a SmartRedis dataset into the database
  procedure :: put_dataset
  !> Retrieve a SmartRedis dataset from the database
  procedure :: get_dataset
  !> Rename the dataset within the database
  procedure :: rename_dataset
  !> Copy a dataset stored in the database into another key
  procedure :: copy_dataset
  !> Delete the dataset from the database
  procedure :: delete_dataset

  procedure :: use_tensor_ensemble_prefix
  procedure :: use_model_ensemble_prefix
  procedure :: set_data_source


  ! Private procedures
  procedure, private :: put_tensor_i8
  procedure, private :: put_tensor_i16
  procedure, private :: put_tensor_i32
  procedure, private :: put_tensor_i64
  procedure, private :: put_tensor_float
  procedure, private :: put_tensor_double
  procedure, private :: unpack_tensor_i8
  procedure, private :: unpack_tensor_i16
  procedure, private :: unpack_tensor_i32
  procedure, private :: unpack_tensor_i64
  procedure, private :: unpack_tensor_float
  procedure, private :: unpack_tensor_double

end type client_type

type dataset_type

end type dataset_type

contains

!> Initializes a new instance of a SmartRedis client
subroutine initialize( this, cluster )
  class(client_type) :: this
  logical, optional :: cluster !< If true, client uses a database cluster (Default: .false.)

end subroutine initialize

logical function isinitialized(this)
  class(client_type) :: this
  isinitialized = .true.
end function isinitialized

!> A destructor for the SmartRedis client
subroutine destructor( this )
  class(client_type) :: this

end subroutine destructor

!> Check if the specified key exists in the database
logical function key_exists(this, key)
  class(client_type)    :: this
  character(len=*)      :: key

end function key_exists

!> Check if the specified model exists in the database
logical function model_exists(this, model_name)
  class(client_type)    :: this
  character(len=*)      :: model_name

end function model_exists

!> Check if the specified tensor exists in the database
logical function tensor_exists(this, tensor_name)
  class(client_type)    :: this
  character(len=*)      :: tensor_name

end function tensor_exists

!> Repeatedly poll the database until the tensor exists or the number of tries is exceeded
logical function poll_tensor( this, tensor_name, poll_frequency_ms, num_tries )
  class(client_type)    :: this
  character(len=*) :: tensor_name !< Key in the database to poll
  integer          :: poll_frequency_ms !< Frequency at which to poll the database (ms)
  integer          :: num_tries !< Number of times to poll the database before failing

end function poll_tensor

!> Repeatedly poll the database until the model exists or the number of tries is exceeded
logical function poll_model( this, model_name, poll_frequency_ms, num_tries )
  class(client_type)    :: this
  character(len=*) :: model_name !< Key in the database to poll
  integer          :: poll_frequency_ms !< Frequency at which to poll the database (ms)
  integer          :: num_tries !< Number of times to poll the database before failing

end function poll_model

!> Repeatedly poll the database until the key exists or the number of tries is exceeded
logical function poll_key( this, key, poll_frequency_ms, num_tries )
  class(client_type)    :: this
  character(len=*) :: key !< Key in the database to poll
  integer          :: poll_frequency_ms !< Frequency at which to poll the database (ms)
  integer          :: num_tries !< Number of times to poll the database before failing

end function poll_key

!> Put a tensor whose Fortran type is the equivalent 'int8' C-type
subroutine put_tensor_i8(this, key, data, dims)
  integer(kind=c_int8_t), dimension(*), target, intent(in) :: data !< Data to be sent
  include 'put_tensor_methods_common.inc'
end subroutine put_tensor_i8

!> Put a tensor whose Fortran type is the equivalent 'int16' C-type
subroutine put_tensor_i16(this, key, data, dims)
  integer(kind=c_int16_t), dimension(*), target, intent(in) :: data !< Data to be sent
  include 'put_tensor_methods_common.inc'
end subroutine put_tensor_i16

!> Put a tensor whose Fortran type is the equivalent 'int32' C-type
subroutine put_tensor_i32(this, key, data, dims)
  integer(kind=c_int32_t), dimension(*), target, intent(in) :: data !< Data to be sent
  include 'put_tensor_methods_common.inc'
end subroutine put_tensor_i32

!> Put a tensor whose Fortran type is the equivalent 'int64' C-type
subroutine put_tensor_i64(this, key, data, dims)
  integer(kind=c_int64_t), dimension(*), target, intent(in) :: data !< Data to be sent
  include 'put_tensor_methods_common.inc'
end subroutine put_tensor_i64

!> Put a tensor whose Fortran type is the equivalent 'float' C-type
subroutine put_tensor_float(this, key, data, dims)
  real(kind=c_float), dimension(*), target, intent(in) :: data !< Data to be sent
  include 'put_tensor_methods_common.inc'
end subroutine put_tensor_float

!> Put a tensor whose Fortran type is the equivalent 'double' C-type
subroutine put_tensor_double(this, key, data, dims)
  real(kind=c_double), dimension(*), target, intent(in) :: data !< Data to be sent
  include 'put_tensor_methods_common.inc'
end subroutine put_tensor_double

!> Put a tensor whose Fortran type is the equivalent 'int8' C-type
subroutine unpack_tensor_i8(this, key, result, dims)
  integer(kind=c_int8_t), dimension(*), target, intent(out) :: result !< Data to be sent
  include 'unpack_tensor_methods_common.inc'
  ! Define the type and call the C-interface
end subroutine unpack_tensor_i8

!> Put a tensor whose Fortran type is the equivalent 'int16' C-type
subroutine unpack_tensor_i16(this, key, result, dims)
  integer(kind=c_int16_t), dimension(*), target, intent(out) :: result !< Data to be sent
  include 'unpack_tensor_methods_common.inc'
end subroutine unpack_tensor_i16

!> Put a tensor whose Fortran type is the equivalent 'int32' C-type
subroutine unpack_tensor_i32(this, key, result, dims)
  integer(kind=c_int32_t), dimension(*), target, intent(out) :: result !< Data to be sent
  include 'unpack_tensor_methods_common.inc'
end subroutine unpack_tensor_i32

!> Put a tensor whose Fortran type is the equivalent 'int64' C-type
subroutine unpack_tensor_i64(this, key, result, dims)
  integer(kind=c_int64_t), dimension(*), target, intent(out) :: result !< Data to be sent
  include 'unpack_tensor_methods_common.inc'
end subroutine unpack_tensor_i64

!> Put a tensor whose Fortran type is the equivalent 'float' C-type
subroutine unpack_tensor_float(this, key, result, dims)
  real(kind=c_float), dimension(*), target, intent(out) :: result !< Data to be sent
  include 'unpack_tensor_methods_common.inc'
end subroutine unpack_tensor_float

!> Put a tensor whose Fortran type is the equivalent 'double' C-type
subroutine unpack_tensor_double(this, key, result, dims)
  real(kind=c_double), dimension(*), target, intent(out) :: result !< Data to be sent
  include 'unpack_tensor_methods_common.inc'
end subroutine unpack_tensor_double

!> Move a tensor to a new key
subroutine rename_tensor(this, key, new_key)
  class(client_type), intent(in) :: this    !< The initialized Fortran SmartRedis client
  character(len=*), intent(in) :: key     !< The key to use to place the tensor
                                          !! excluding null terminating character
  character(len=*), intent(in) :: new_key !< The new tensor key

end subroutine rename_tensor

!> Delete a tensor
subroutine delete_tensor(this, key)
  class(client_type), intent(in) :: this !<  The initialized Fortran SmartRedis client
  character(len=*), intent(in) :: key  !< The key to use to place the tensor
end subroutine delete_tensor

!> Copy a tensor to the destination key
subroutine copy_tensor(this, src_name, dest_name)
  class(client_type), intent(in) :: this      !< The initialized Fortran SmartRedis client
  character(len=*), intent(in) :: src_name  !< The key to use to place the tensor
                                            !! excluding null terminating character
  character(len=*), intent(in) :: dest_name !< The new tensor key

end subroutine copy_tensor

!> Retrieve the model from the database
subroutine get_model(this, key, model)
  class(client_type),               intent(in   ) :: this  !< An initialized SmartRedis client
  character(len=*),                 intent(in   ) :: key   !< The key to use to place the model
  character(len=*),                 intent(  out) :: model !< The model as a continuous buffer

end subroutine get_model

!> Load the machine learning model from a file and set the configuration
subroutine set_model_from_file( this, key, model_file, backend, device, batch_size, min_batch_size, tag, &
    inputs, outputs )
  class(client_type),             intent(in) :: this           !< An initialized SmartRedis client
  character(len=*),               intent(in) :: key            !< The key to use to place the model
  character(len=*),               intent(in) :: model_file     !< The file storing the model
  character(len=*),               intent(in) :: backend        !< The name of the backend (TF, TFLITE, TORCH, ONNX)
  character(len=*),               intent(in) :: device         !< The name of the device (CPU, GPU, GPU:0, GPU:1...)
  integer,                        optional,  intent(in) :: batch_size     !< The batch size for model execution
  integer,                        optional,  intent(in) :: min_batch_size !< The minimum batch size for model execution
  character(len=*),               optional,  intent(in) :: tag            !< A tag to attach to the model for
                                                                          !! information purposes
  character(len=*), dimension(:), optional,  intent(in) :: inputs         !< One or more names of model input nodes (TF
                                                                          !! models)
  character(len=*), dimension(:), optional,  intent(in) :: outputs        !< One or more names of model output nodes (TF models)

end subroutine set_model_from_file

subroutine set_model( this, key, model, backend, device, batch_size, min_batch_size, tag, &
    inputs, outputs )
  class(client_type),             intent(in) :: this           !< An initialized SmartRedis client
  character(len=*),               intent(in) :: key            !< The key to use to place the model
  character(len=*),               intent(in) :: model          !< The binary representaiton o
  character(len=*),               intent(in) :: backend        !< The name of the backend (TF, TFLITE, TORCH, ONNX)
  character(len=*),               intent(in) :: device         !< The name of the device (CPU, GPU, GPU:0, GPU:1...)
  integer,                        intent(in) :: batch_size     !< The batch size for model execution
  integer,                        intent(in) :: min_batch_size !< The minimum batch size for model execution
  character(len=*),               intent(in) :: tag            !< A tag to attach to the model for information purposes
  character(len=*), dimension(:), intent(in) :: inputs         !< One or more names of model input nodes (TF models)
  character(len=*), dimension(:), intent(in) :: outputs        !< One or more names of model output nodes (TF models)

end subroutine set_model

subroutine run_model(this, key, inputs, outputs)
  class(client_type),             intent(in) :: this           !< An initialized SmartRedis client
  character(len=*),               intent(in) :: key            !< The key to use to place the model
  character(len=*), dimension(:), intent(in) :: inputs         !< One or more names of model input nodes (TF models)
  character(len=*), dimension(:), intent(in) :: outputs        !< One or more names of model output nodes (TF models)

end subroutine run_model

!> Retrieve the script from the database
subroutine get_script(this, key, script)
  class(client_type),               intent(in   ) :: this  !< An initialized SmartRedis client
  character(len=*),                 intent(in   ) :: key   !< The key to use to place the script
  character(len=*),                 intent(  out) :: script !< The script as a continuous buffer

end subroutine get_script

subroutine set_script_from_file( this, key, device, script_file )
  class(client_type),             intent(in) :: this        !< An initialized SmartRedis client
  character(len=*),               intent(in) :: key         !< The key to use to place the script
  character(len=*),               intent(in) :: device      !< The name of the device (CPU, GPU, GPU:0, GPU:1...)
  character(len=*),               intent(in) :: script_file !< The file storing the script

end subroutine set_script_from_file

subroutine set_script( this, key, device, script )
  class(client_type),             intent(in) :: this   !< An initialized SmartRedis client
  character(len=*),               intent(in) :: key    !< The key to use to place the script
  character(len=*),               intent(in) :: device !< The name of the device (CPU, GPU, GPU:0, GPU:1...)
  character(len=*),               intent(in) :: script !< The file storing the script

end subroutine set_script

subroutine run_script(this, key, func, inputs, outputs)
  class(client_type),             intent(in) :: this           !< An initialized SmartRedis client
  character(len=*),               intent(in) :: key            !< The key to use to place the script
  character(len=*),               intent(in) :: func           !< The name of the function in the script to call
  character(len=*), dimension(:), intent(in) :: inputs         !< One or more names of script input nodes (TF scripts)
  character(len=*), dimension(:), intent(in) :: outputs        !< One or more names of script output nodes (TF scripts)

end subroutine run_script
!> Store a dataset in the database
subroutine put_dataset( this, dataset )
  class(client_type), intent(in) :: this
  type(dataset_type), intent(in) :: dataset !< Dataset to store in the dataset

end subroutine put_dataset

!> Retrieve a dataset from the database
type(dataset_type) function get_dataset( this, name )
  class(client_type) :: this
  character(len=*)   :: name !< Name of the dataset to get

end function get_dataset

!> Rename a dataset stored in the database
subroutine rename_dataset( this, name, new_name )
  class(client_type) :: this
  character(len=*)   :: name     !< Original name of the dataset
  character(len=*)   :: new_name !< New name of the dataset

end subroutine rename_dataset

!> Copy a dataset within the database to a new name
subroutine copy_dataset( this, name, new_name )
  class(client_type) :: this
  character(len=*)   :: name     !< Original name of the dataset
  character(len=*)   :: new_name !< New name of the dataset

end subroutine copy_dataset

!> Delete a dataset stored within a database
subroutine delete_dataset( this, name )
  class(client_type) :: this
  character(len=*)   :: name !< Name of the dataset to delete

end subroutine delete_dataset

!> Set the data source (i.e. key prefix for get functions)
subroutine set_data_source( this, source_id )
  class(client_type), intent(in) :: this
  character(len=*),   intent(in) :: source_id


end subroutine set_data_source

!> Set whether names of model and script entities should be prefixed (e.g. in an ensemble) to form database keys.
!! Prefixes will only be used if they were previously set through the environment variables SSKEYOUT and SSKEYIN.
!! Keys of entities created before this function is called will not be affected. By default, the client does not
!! prefix model and script keys.
subroutine use_model_ensemble_prefix( this, use_prefix )
  class(client_type),   intent(in) :: this
  logical,              intent(in) :: use_prefix

end subroutine use_model_ensemble_prefix


!> Set whether names of tensor and dataset entities should be prefixed (e.g. in an ensemble) to form database keys.
!! Prefixes will only be used if they were previously set through the environment variables SSKEYOUT and SSKEYIN.
!! Keys of entities created before this function is called will not be affected. By default, the client prefixes
!! tensor and dataset keys with the first prefix specified with the SSKEYIN and SSKEYOUT environment variables.
subroutine use_tensor_ensemble_prefix( this, use_prefix )
  class(client_type),   intent(in) :: this
  logical,              intent(in) :: use_prefix

end subroutine use_tensor_ensemble_prefix

end module smartredis_client
