module esmflds
  use ESMF, only                  : ESMF_SUCCESS, ESMF_FAILURE, ESMF_LOGWRITE
  use med_kind_mod, only          : CX=>SHR_KIND_CX, CS=>SHR_KIND_CS, CL=>SHR_KIND_CL, R8=>SHR_KIND_R8
  use med_internalstate_mod, only : compname, compocn, compatm, compice, comprof
  use med_internalstate_mod, only : mapfcopy, mapnames, mapunset
  use med_utils_mod        , only : chkerr => med_utils_ChkErr
  use shr_log_mod          , only : shr_log_error
  implicit none
  private

  !-----------------------------------------------
  ! PUblic methods
  !-----------------------------------------------

  public :: med_fldList_init1

  public :: med_fldList_addfld_from
  public :: med_fldList_addmap_from
  public :: med_fldList_addfld_to
  public :: med_fldList_addmrg_to

  public :: med_fldList_addfld_ocnalb
  public :: med_fldList_addmap_ocnalb

  public :: med_fldList_addfld_aoflux
  public :: med_fldList_addmap_aoflux

  private :: med_fldList_AddFld
  private :: med_fldList_AddMap
  private :: med_fldList_AddMrg
  public :: med_fldList_findName
  public :: med_fldList_GetFldNames
  public :: med_fldList_GetNumFlds
  public :: med_fldList_GetFldInfo
  public :: med_fld_GetFldInfo
  public :: med_fldList_Realize
  public :: med_fldList_Document_Mapping
  public :: med_fldList_Document_Merging
  public :: med_fldList_GetFldListFr
  public :: med_fldList_GetFldListTo
  public :: med_fldList_GetaofluxFldList
  public :: med_fldList_GetocnalbFldList
  !-----------------------------------------------
  ! Types and instantiations that determine fields, mappings, mergings
  !-----------------------------------------------

  type, public :: med_fldList_entry_type
     character(CS) :: stdname
     character(CS) :: shortname

     ! Mapping fldsFr data - for mediator import fields
     integer      , allocatable :: mapindex(:)
     character(CS), allocatable :: mapnorm(:)
     character(CX), allocatable :: mapfile(:)

     ! Merging fldsTo data - for mediator export field
     character(CS), allocatable :: merge_fields(:)
     character(CS), allocatable :: merge_types(:)
     character(CS), allocatable :: merge_fracnames(:)
     type(med_fldList_entry_type), pointer :: next => null()
  end type med_fldList_entry_type

  ! The above would be the field name to merge from
  ! e.g. for Sa_z in lnd
  !    merge_field(compatm) = 'Sa_z'
  !    merge_type(comptm) = 'copy'  (could also have 'copy_with_weighting')

  type, public :: med_fldList_type
     type (med_fldList_entry_type) :: fields
  end type med_fldList_type

  !-----------------------------------------------
  ! Instantiate derived types
  !-----------------------------------------------
  type (med_fldList_type), allocatable, target :: fldListTo(:) ! advertise fields to components
  type (med_fldList_type), allocatable, target :: fldListFr(:) ! advertise fields from components

  type (med_fldList_type), target :: fldlist_aoflux
  type (med_fldList_type), target :: fldlist_ocnalb

  integer                    :: rc
  character(len=CL)          :: infostr
  character(len=*),parameter :: u_FILE_u = &
     __FILE__

!================================================================================
contains
!================================================================================

  subroutine med_fldlist_init1(ncomps)
    integer, intent(in) :: ncomps
    allocate(fldlistTo(ncomps))
    allocate(fldlistFr(ncomps))
  end subroutine med_fldlist_init1

  !================================================================================

  function med_fldList_GetaofluxFldList() result(fldList)
    ! Return a pointer to the aoflux fldlist
    type(med_fldList_type), pointer :: fldList

    fldList => fldlist_aoflux
  end function Med_FldList_GetaofluxFldList

  !================================================================================

  function med_fldList_GetocnalbFldList() result(fldList)
    ! Return a pointer to the ocnalb fldlist
    type(med_fldList_type), pointer :: fldList

    fldList => fldlist_ocnalb
  end function Med_FldList_GetocnalbFldList

  !================================================================================

  function med_fldList_GetFldListFr(index) result(fldList)
    ! Return a pointer to the FldListFr(index)
    integer, intent(in) :: index
    type(med_fldList_type), pointer :: fldList

    fldList => fldListFr(index)
  end function Med_FldList_GetFldListFr

  !================================================================================

  function med_fldList_GetFldListTo(index) result(fldList)
    ! Return a pointer to the FldListTo(index)
    integer, intent(in) :: index
    type(med_fldList_type), pointer :: fldList

    fldList => fldListTo(index)
  end function Med_FldList_GetFldListTo

  !================================================================================

  subroutine med_fldList_addfld_from(index, stdname, shortname)
    ! add a fld with name stdname to the FldListFr list
    integer, intent(in) :: index
    character(len=*)             , intent(in)             :: stdname
    character(len=*)             , intent(in)  , optional :: shortname

    call med_fldList_AddFld(FldListFr(index)%fields, stdname, shortname)

  end subroutine med_fldList_addfld_from

  !================================================================================

  subroutine med_fldList_addfld_aoflux(stdname, shortname)
    ! add a fld to the aoflux fldList
    character(len=*)             , intent(in)             :: stdname
    character(len=*)             , intent(in)  , optional :: shortname

    call med_fldList_AddFld(fldlist_aoflux%fields, stdname, shortname)

  end subroutine med_fldList_addfld_aoflux

  !================================================================================

  subroutine med_fldList_addfld_ocnalb(stdname, shortname)
    character(len=*)             , intent(in)             :: stdname
    character(len=*)             , intent(in)  , optional :: shortname

    call med_fldList_AddFld(fldlist_ocnalb%fields, stdname, shortname)

  end subroutine med_fldList_addfld_ocnalb

  !================================================================================

  subroutine med_fldList_addfld_to(index, stdname, shortname)
    integer, intent(in) :: index
    character(len=*)             , intent(in)             :: stdname
    character(len=*)             , intent(in)  , optional :: shortname

    call med_fldList_AddFld(FldListTo(index)%fields, stdname, shortname)

  end subroutine med_fldList_addfld_to

  !================================================================================

  subroutine med_fldList_findName(fields, stdname, found, lastfld)
    ! on return if found == .true. lastfield is the field matching stdname
    ! if found == .false. lastfield is the last field in the list
    type(med_fldList_entry_type) , intent(in), target           :: fields
    character(len=*)             , intent(in)             :: stdname
    logical                      , intent(out)            :: found
    type(med_fldList_entry_type) , intent(out), pointer                :: lastfld

    lastfld => fields
    found = .false.
    do while(associated(lastfld%next))
       if (trim(stdname) == trim(lastfld%stdname)) exit
       lastfld => lastfld%next
    enddo
    ! Check the lastfld
    if (trim(stdname) == trim(lastfld%stdname)) found = .true.

  end subroutine med_fldList_findName

  !================================================================================

  subroutine med_fldList_AddFld(fields, stdname, shortname)
    ! ----------------------------------------------
    ! Add an entry to to the flds array
    ! Use pointers to create an extensible allocatable array.
    ! to allow the size of flds to grow, the process for
    ! adding a new field is:
    ! 1) allocate newflds to be N (one element larger than flds)
    ! 2) copy flds into first N-1 elements of newflds
    ! 3) newest flds entry is Nth element of newflds
    ! 4) deallocate / nullify flds
    ! 5) point flds => newflds
    ! ----------------------------------------------

    type(med_fldList_entry_type) , target                :: fields
    character(len=*)             , intent(in)             :: stdname
    character(len=*)             , intent(in)  , optional :: shortname

    ! local variables
    logical :: found
    integer :: mapsize, mrgsize
    type(med_fldList_entry_type), pointer :: newfld
    character(len=*), parameter :: subname='(med_fldList_AddFld)'
    ! ----------------------------------------------

    call med_fldList_findName(fields, stdname, found, newfld)
    ! create new entry if fldname is not in original list
    mapsize = size(fldListTo)
    mrgsize = size(fldListFr)

    if (.not. found) then
       ! 1) allocate newfld to be size (one element larger than input flds)
       ! the if statement allows the first entry to be filed
       if(allocated(newfld%mapindex)) then
          allocate(newfld%next)
          newfld => newfld%next
       endif

       ! 2) now update flds information for new entry
       newfld%stdname   = trim(stdname)
       if (present(shortname)) then
          newfld%shortname = trim(shortname)
       else
          newfld%shortname = trim(stdname)
       end if
       allocate(newfld%mapindex(mapsize))
       allocate(newfld%mapnorm(mapsize))
       allocate(newfld%mapfile(mapsize))
       allocate(newfld%merge_fields(mrgsize))
       allocate(newfld%merge_types(mrgsize))
       allocate(newfld%merge_fracnames(mrgsize))
       newfld%mapindex(:) = mapunset
       newfld%mapnorm(:) = 'unset'
       newfld%mapfile(:) = 'unset'
       newfld%merge_fields(:) = 'unset'
       newfld%merge_types(:) = 'unset'
       newfld%merge_fracnames(:) = 'unset'
    end if

  end subroutine med_fldList_AddFld

  !================================================================================

  subroutine med_fldList_addmrg_to(index, fldname, mrg_from, mrg_fld, mrg_type, mrg_fracname, rc)

    ! ----------------------------------------------
    ! Determine mrg entry or entries in flds aray
    ! ----------------------------------------------

    ! input/output variables
    integer                      , intent(in)           :: index
    character(len=*)             , intent(in)           :: fldname
    integer                      , intent(in)           :: mrg_from
    character(len=*)             , intent(in)           :: mrg_fld
    character(len=*)             , intent(in)           :: mrg_type
    character(len=*)             , intent(in), optional :: mrg_fracname
    integer                      , intent(out), optional :: rc

    call med_FldList_addMrg(fldListTo(index)%fields, fldname, mrg_from, mrg_fld, mrg_type, mrg_fracname)

  end subroutine med_fldList_addmrg_to

  !================================================================================

  subroutine med_fldList_AddMrg(flds, fldname, mrg_from, mrg_fld, mrg_type, mrg_fracname)

    ! ----------------------------------------------
    ! Determine mrg entry or entries in flds aray
    ! ----------------------------------------------

    ! input/output variables
    type(med_fldList_entry_type) , intent(in), target   :: flds
    character(len=*)             , intent(in)           :: fldname
    integer                      , intent(in)           :: mrg_from
    character(len=*)             , intent(in)           :: mrg_fld
    character(len=*)             , intent(in)           :: mrg_type
    character(len=*)             , intent(in), optional :: mrg_fracname

    ! local variables
    integer :: rc
    type(med_fldList_entry_type), pointer :: newfld
    character(len=*), parameter :: subname='(med_fldList_AddMrg)'
    ! ----------------------------------------------

    newfld => med_fldList_GetFld(flds, fldname, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    newfld%merge_fields(mrg_from) = mrg_fld
    newfld%merge_types(mrg_from) = mrg_type
    if (present(mrg_fracname)) then
       newfld%merge_fracnames(mrg_from) = mrg_fracname
    end if

  end subroutine med_fldList_AddMrg

  !================================================================================

  function med_fldList_GetFld(fields, fldname, rc) result(newfld)
    use ESMF, only : ESMF_LogWrite, ESMF_LOGMSG_INFO

    type(med_fldList_entry_type) , intent(in), target :: fields
    character(len=*)                  , intent(in)    :: fldname

    type(med_fldList_entry_type), pointer :: newfld
    logical :: found
    integer :: rc
    character(len=*), parameter :: subname='(med_fldList_GetFld)'


    call med_fldList_findName(fields, fldname, found, newfld)

    rc = ESMF_SUCCESS
    if(.not. found) then
       rc = ESMF_FAILURE
       newfld => fields
       do while(associated(newfld))
          write(6,*) trim(subname)//' input flds entry is ',trim(newfld%stdname)
          newfld => newfld%next
       end do
       call shr_log_error(subname // 'ERROR: fldname '// trim(fldname) // ' not found in input flds', rc=rc)
       return
    endif

  end function med_fldList_GetFld

  !================================================================================

  subroutine med_fldList_addmap_from(index, fldname, destcomp, maptype, mapnorm, mapfile)
    integer, intent(in) :: index
    character(len=*)                  , intent(in)    :: fldname
    integer                            , intent(in)    :: destcomp
    integer                            , intent(in)    :: maptype
    character(len=*)                   , intent(in)    :: mapnorm
    character(len=*), optional         , intent(in)    :: mapfile

    call med_fldList_AddMap(FldListFr(index)%fields, fldname, destcomp, maptype, mapnorm, mapfile)

  end subroutine med_fldList_addmap_from

  !================================================================================

  subroutine med_fldList_addmap_aoflux(fldname, destcomp, maptype, mapnorm, mapfile)
    character(len=*)                  , intent(in)    :: fldname
    integer                            , intent(in)    :: destcomp
    integer                            , intent(in)    :: maptype
    character(len=*)                   , intent(in)    :: mapnorm
    character(len=*), optional         , intent(in)    :: mapfile

    call med_fldList_AddMap(fldlist_aoflux%fields, fldname, destcomp, maptype, mapnorm, mapfile)

  end subroutine med_fldList_addmap_aoflux

  !================================================================================

  subroutine med_fldList_addmap_ocnalb(fldname, destcomp, maptype, mapnorm, mapfile)
    character(len=*)                  , intent(in)    :: fldname
    integer                            , intent(in)    :: destcomp
    integer                            , intent(in)    :: maptype
    character(len=*)                   , intent(in)    :: mapnorm
    character(len=*), optional         , intent(in)    :: mapfile

    call med_fldList_AddMap(fldlist_ocnalb%fields, fldname, destcomp, maptype, mapnorm, mapfile)

  end subroutine med_fldList_addmap_ocnalb

  !================================================================================

  subroutine med_fldList_AddMap(fields, fldname, destcomp, maptype, mapnorm, mapfile)
    ! intput/output variables
    type(med_fldList_entry_type) , intent(in), target :: fields
    character(len=*)                  , intent(in)    :: fldname
    integer                            , intent(in)    :: destcomp
    integer                            , intent(in)    :: maptype
    character(len=*)                   , intent(in)    :: mapnorm
    character(len=*), optional         , intent(in)    :: mapfile

    ! local variables
    type(med_fldList_entry_type), pointer :: newfld
    integer :: rc

    character(len=CX)                                  :: lmapfile
    character(len=*),parameter  :: subname='(med_fldList_AddMap)'
    ! ----------------------------------------------
    lmapfile = 'unset'
    if (present(mapfile)) lmapfile = mapfile

    newfld => med_fldList_GetFld(fields, fldname, rc)
    if (chkerr(rc,__LINE__,u_FILE_u)) return
    ! Note - default values are already set for the fld entries - so only non-default
    ! values need to be set below
    ! If mapindex is mapfcopy - create a redistribution route handle
    ! If mapfile is idmap - create a redistribution route nhandle
    ! If mapfile is unset then create the mapping route handle at run time

    newfld%mapindex(destcomp) = maptype
    newfld%mapnorm(destcomp)  = trim(mapnorm)
    newfld%mapfile(destcomp)  = trim(lmapfile)

    ! overwrite values if appropriate
    if (newfld%mapindex(destcomp) == mapfcopy) then
       newfld%mapfile(destcomp) = 'unset'
       newfld%mapnorm(destcomp) = 'unset'
    else if (trim(newfld%mapfile(destcomp)) == 'idmap') then
       newfld%mapindex(destcomp) = mapfcopy
       newfld%mapnorm(destcomp) = 'unset'
    end if

  end subroutine med_fldList_AddMap

  !================================================================================

  subroutine med_fldList_Realize(state, fldList, flds_scalar_name, flds_scalar_num, &
       grid, mesh, tag, rc)

    use NUOPC             , only : NUOPC_GetStateMemberLists, NUOPC_IsConnected, NUOPC_Realize
    use NUOPC             , only : NUOPC_GetAttribute
    use ESMF              , only : ESMF_MeshLoc_Element, ESMF_FieldCreate, ESMF_TYPEKIND_R8
    use ESMF              , only : ESMF_MAXSTR, ESMF_Field, ESMF_State, ESMF_Grid, ESMF_Mesh
    use ESMF              , only : ESMF_StateGet, ESMF_LogFoundError
    use ESMF              , only : ESMF_LogWrite, ESMF_FAILURE, ESMF_LOGERR_PASSTHRU
    use ESMF              , only : ESMF_LOGMSG_INFO, ESMF_StateRemove, ESMF_SUCCESS
    use ESMF              , only : ESMF_STATEINTENT_IMPORT, ESMF_STATEINTENT_EXPORT, ESMF_StateIntent_Flag
    use ESMF              , only : ESMF_RC_ARG_BAD, ESMF_LogSetError, operator(==)
    ! input/output variables
    type(ESMF_State)            , intent(inout)            :: state
    type(med_fldlist_type)      , intent(in), target       :: fldList
    character(len=*)            , intent(in)               :: flds_scalar_name
    integer                     , intent(in)               :: flds_scalar_num
    character(len=*)            , intent(in)               :: tag
    integer                     , intent(inout)            :: rc
    type(ESMF_Grid)             , intent(in)    , optional :: grid
    type(ESMF_Mesh)             , intent(in)    , optional :: mesh

    ! local variables
    type(med_fldList_entry_type), pointer :: newfld
    integer                         :: itemCount
    integer                         :: n
    type(ESMF_Field)                :: field
    character(CS)                   :: shortname
    character(ESMF_MAXSTR)          :: transferActionAttr
    type(ESMF_StateIntent_Flag)     :: stateIntent
    character(ESMF_MAXSTR)          :: transferAction
    character(ESMF_MAXSTR), pointer :: StandardNameList(:)
    character(ESMF_MAXSTR), pointer :: ConnectedList(:)
    character(ESMF_MAXSTR), pointer :: NameSpaceList(:)
    character(ESMF_MAXSTR), pointer :: itemNameList(:)
    character(len=*),parameter  :: subname='(med_fldList_Realize)'
    ! ----------------------------------------------

    rc = ESMF_SUCCESS

    if (present(grid) .and. present(mesh)) then
       call shr_log_error(trim(subname)//trim(tag)//": ERROR both grid and mesh not allowed", &
            line=__LINE__, file=u_FILE_u, rc=rc)
       return
    endif

    nullify(StandardNameList)
    nullify(ConnectedList)
    nullify(NameSpaceList)
    nullify(ItemNameList)

    call ESMF_StateGet(state, itemCount=itemCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    write(infostr,'(i6)') itemCount
    call ESMF_LogWrite(trim(subname)//trim(tag)//" count = "//trim(infostr), ESMF_LOGMSG_INFO)
    if (itemCount > 0) then
       allocate(itemNameList(itemCount))
       call ESMF_StateGet(state, itemNameList=itemNameList, rc=rc)
       if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
       do n = 1,itemCount
          call ESMF_LogWrite(trim(subname)//trim(tag)//" itemNameList = "//trim(itemNameList(n)), ESMF_LOGMSG_INFO)
       enddo
       deallocate(itemNameList)
    endif

#if (1 == 0)
    call NUOPC_GetStateMemberLists(state, StandardNameList=StandardNameList, ConnectedList=ConnectedList, &
         NamespaceList=NamespaceList, itemNameList=itemNameList, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
    write(infostr,'(i6)') size(StandardNameList)
    call ESMF_LogWrite(trim(subname)//trim(tag)//" size = "//trim(infostr), ESMF_LOGMSG_INFO)

    do n = 1,size(StandardNameList)
       call ESMF_LogWrite(trim(subname)//trim(tag)//" StandardNameList = "//trim(StandardNameList(n)), &
            ESMF_LOGMSG_INFO)
    enddo
    do n = 1,size(ConnectedList)
       call ESMF_LogWrite(trim(subname)//trim(tag)//" ConnectedList = "//trim(ConnectedList(n)), &
            ESMF_LOGMSG_INFO)
    enddo
    do n = 1,size(NamespaceList)
       call ESMF_LogWrite(trim(subname)//trim(tag)//" NamespaceList = "//trim(NamespaceList(n)), &
            ESMF_LOGMSG_INFO)
    enddo
    do n = 1,size(ItemnameList)
       call ESMF_LogWrite(trim(subname)//trim(tag)//" ItemnameList = "//trim(ItemnameList(n)), &
            ESMF_LOGMSG_INFO)
    enddo
#endif

    call ESMF_StateGet(state, stateIntent=stateIntent, rc=rc)
    if (stateIntent==ESMF_STATEINTENT_EXPORT) then
       transferActionAttr="ProducerTransferAction"
    elseif (stateIntent==ESMF_STATEINTENT_IMPORT) then
       transferActionAttr="ConsumerTransferAction"
    else
       call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
            msg="The stateIntent must either be IMPORT or EXPORT here.", &
            line=__LINE__, &
            file=__FILE__, &
            rcToReturn=rc)
       return  ! bail out
    endif

    newfld => fldList%fields
    do while(associated(newfld))
       shortname = newfld%shortname

       ! call ESMF_LogWrite(subname//' fld = '//trim(shortname), ESMF_LOGMSG_INFO)
       if (NUOPC_IsConnected(state, fieldName=shortname)) then

          call ESMF_StateGet(state, field=field, itemName=trim(shortname), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

          call NUOPC_GetAttribute(field, name=TransferActionAttr, value=transferAction, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

          if (trim(transferAction) == "accept") then  ! accept

             call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(shortname)//" is connected, grid/mesh TBD", &
                  ESMF_LOGMSG_INFO)

          else   ! provide

             if (shortname == trim(flds_scalar_name)) then
                call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(shortname)//" is connected on root pe", &
                     ESMF_LOGMSG_INFO)
                call SetScalarField(field, flds_scalar_name, flds_scalar_num, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
             elseif (present(grid)) then
                call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(shortname)//" is connected using grid", &
                     ESMF_LOGMSG_INFO)
                ! Create the field
                field = ESMF_FieldCreate(grid, ESMF_TYPEKIND_R8, name=shortname,rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
             elseif (present(mesh)) then
                call ESMF_LogWrite(trim(subname)//trim(tag)//" Field = "//trim(shortname)//" is connected using mesh", &
                     ESMF_LOGMSG_INFO)
                ! Create the field
                field = ESMF_FieldCreate(mesh, ESMF_TYPEKIND_R8, name=shortname, meshloc=ESMF_MESHLOC_ELEMENT, rc=rc)
                if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return
             else
                call shr_log_error(trim(subname)//trim(tag)//": ERROR grid or mesh expected", &
                     line=__LINE__, file=u_FILE_u, rc=rc)
                return
             endif

             ! NOW call NUOPC_Realize
             call NUOPC_Realize(state, field=field, rc=rc)
             if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

          endif

       else

          call ESMF_LogWrite(subname // trim(tag) // " Field = "// trim(shortname) // " is not connected.", &
               ESMF_LOGMSG_INFO)
          call ESMF_StateRemove(state, (/shortname/), rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

       end if
       newfld => newfld%next
    end do

    call ESMF_LogWrite(subname//' done ', ESMF_LOGMSG_INFO)

  contains  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    subroutine SetScalarField(field, flds_scalar_name, flds_scalar_num, rc)
      ! ----------------------------------------------
      ! create a field with scalar data on the root pe
      ! ----------------------------------------------
      use ESMF, only : ESMF_Field, ESMF_DistGrid, ESMF_Grid
      use ESMF, only : ESMF_DistGridCreate, ESMF_GridCreate, ESMF_LogFoundError, ESMF_LOGERR_PASSTHRU
      use ESMF, only : ESMF_FieldCreate, ESMF_GridCreate, ESMF_TYPEKIND_R8
      type(ESMF_Field) , intent(inout) :: field
      character(len=*) , intent(in)    :: flds_scalar_name
      integer          , intent(in)    :: flds_scalar_num
      integer          , intent(inout) :: rc

      ! local variables
      type(ESMF_Distgrid) :: distgrid
      type(ESMF_Grid)     :: grid
      character(len=*), parameter :: subname='(SetScalarField)'
      ! ----------------------------------------------

      rc = ESMF_SUCCESS

      ! create a DistGrid with a single index space element, which gets mapped onto DE 0.
      distgrid = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      grid = ESMF_GridCreate(distgrid, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

      field = ESMF_FieldCreate(name=trim(flds_scalar_name), &
           grid=grid, &
           typekind=ESMF_TYPEKIND_R8, &
           ungriddedLBound=(/1/), &
           ungriddedUBound=(/flds_scalar_num/), &
           gridToFieldMap=(/2/), &
           rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, line=__LINE__, file=u_FILE_u)) return

    end subroutine SetScalarField

  end subroutine med_fldList_Realize

  !================================================================================

  subroutine med_fldList_GetFldInfo(fldList, fldindex, compsrc, stdname, shortname, mapindex, mapFile, mapnorm, merge_fields, merge_type, merge_fracname, rc)
    ! ----------------------------------------------
    ! Get field info
    ! ----------------------------------------------
    type(med_fldList_type) , intent(in), target          :: fldList
    integer                      , intent(in)            :: fldindex
    integer                      , optional, intent(in)  :: compsrc
    integer                      , optional, intent(out) :: mapindex
    character(len=*)             , optional, intent(out) :: mapfile
    character(len=*)             , optional, intent(out) :: mapnorm
    character(len=*)             , optional, intent(out) :: stdname
    character(len=*)             , optional, intent(out) :: shortname
    character(len=*)             , optional, intent(out) :: merge_fields
    character(len=*)             , optional, intent(out) :: merge_type
    character(len=*)             , optional, intent(out) :: merge_fracname
    integer                      , optional, intent(out) :: rc

    ! local variables
    type(med_fldList_entry_type), pointer :: newfld
    integer :: i
    integer :: lcompsrc
    character(len=*), parameter :: subname='(med_fldList_GetFldInfo)'
    ! ----------------------------------------------
    i = 0
    lcompsrc = 1
    newfld => fldList%fields
    do while(associated(newfld))
       i = i+1
       if (i==fldindex) exit
       newfld => newfld%next
    enddo
    if( .not. associated(newfld)) then
       call shr_log_error(subname//' No field found', rc=rc)
       return
    endif
    call med_fld_GetFldInfo(newfld, compsrc, stdname, shortname, mapindex, mapFile, mapnorm, merge_fields, merge_type, merge_fracname, rc)

  end subroutine med_fldList_GetFldInfo

  !================================================================================

  subroutine med_fld_GetFldInfo(newfld, compsrc, stdname, shortname, mapindex, mapFile, mapnorm, merge_fields, merge_type, merge_fracname, rc)
    ! ----------------------------------------------
    ! Get field info
    ! ----------------------------------------------
    type(med_fldList_entry_type) , intent(in)            :: newfld
    integer                      , optional, intent(in)  :: compsrc
    integer                      , optional, intent(out) :: mapindex
    character(len=*)             , optional, intent(out) :: mapfile
    character(len=*)             , optional, intent(out) :: mapnorm
    character(len=*)             , optional, intent(out) :: stdname
    character(len=*)             , optional, intent(out) :: shortname
    character(len=*)             , optional, intent(out) :: merge_fields
    character(len=*)             , optional, intent(out) :: merge_type
    character(len=*)             , optional, intent(out) :: merge_fracname
    integer                      , optional, intent(out) :: rc

    ! local variables
    integer :: lrc
    integer :: lcompsrc
    character(len=*), parameter :: subname='(med_fld_GetFldInfo)'
    lrc = ESMF_SUCCESS

    lcompsrc = -1
    if(present(compsrc)) lcompsrc = compsrc

    if(present(stdname)) then
       stdname   = newfld%stdname
    endif
    if(present(shortname)) then
       shortname = newfld%shortname
    endif

    if(present(mapindex)) then
       if(lcompsrc < 0) call shr_log_error("In med_fld_GetFldInfo mapindex requiring compsrc was requested but compsrc was not provided. ", rc=lrc)
       mapindex    = newfld%mapindex(lcompsrc)
    endif
    if(present(mapfile)) then
       if(lcompsrc < 0) call shr_log_error("In med_fld_GetFldInfo mapfile requiring compsrc was requested but compsrc was not provided. ", rc=lrc)
       mapfile    = newfld%mapfile(lcompsrc)
    endif
    if(present(mapnorm)) then
       if(lcompsrc < 0) call shr_log_error("In med_fld_GetFldInfo mapnorm requiring compsrc was requested but compsrc was not provided. ", rc=lrc)
       mapnorm    = newfld%mapnorm(lcompsrc)
    endif
    if(present(merge_fields)) then
       if(lcompsrc < 0) call shr_log_error("In med_fld_GetFldInfo merge_fields requiring compsrc was requested but compsrc was not provided. ", rc=lrc)
       merge_fields    = newfld%merge_fields(lcompsrc)
    endif
    if(present(merge_type)) then
       if(lcompsrc < 0) call shr_log_error("In med_fld_GetFldInfo merge_type requiring compsrc was requested but compsrc was not provided. ", rc=lrc)
       merge_type     = newfld%merge_types(lcompsrc)
    endif
    if(present(merge_fracname)) then
       if(lcompsrc < 0) call shr_log_error("In med_fld_GetFldInfo merge_fracname requiring compsrc was requested but compsrc was not provided. ", rc=lrc)
       merge_fracname = newfld%merge_fracnames(lcompsrc)
    endif
    if(present(rc)) rc=lrc

  end subroutine med_fld_GetFldInfo

  !================================================================================

  integer function med_fldList_GetNumFlds(fldList)

    ! input/output variables
    type(med_fldList_type), intent(in), target  :: fldList
    ! ----------------------------------------------
    type(med_fldList_entry_type), pointer :: newfld

    newfld => fldList%fields
    med_fldList_GetNumFlds = 0
    do while(associated(newfld))
       if(allocated(newfld%mapindex)) then
          med_fldList_GetNumFlds = med_fldList_GetNumFlds + 1
       endif
       newfld => newfld%next
    end do

  end function med_fldList_GetNumFlds

  !================================================================================

  subroutine med_fldList_GetFldNames(fields, fldnames, rc)

    use ESMF, only : ESMF_SUCCESS

    ! input/output variables
    type(med_fldList_entry_type) , intent(in), target     :: fields
    character(len=*)             , intent(inout), pointer   :: fldnames(:)
    integer, optional            , intent(out) :: rc

    !local variables
    type(med_fldList_entry_type), pointer :: newfld
    integer :: n
    character(len=CL) :: msg
    ! ----------------------------------------------

    if(present(rc)) rc = ESMF_SUCCESS
    if (.not. associated(fldnames) .or. .not. allocated(fields%mapindex)) then
       write(msg, *) "med_fldList_GetFldNames: ERROR either fields or fldnames have not been allocated. ",associated(fldnames), allocated(fields%mapindex)
       call shr_log_error(msg, rc=rc)
       return
    endif
    n = 0
    newfld => fields
    do while(associated(newfld))
       n = n+1
       fldnames(n) = trim(newfld%shortname)
       newfld => newfld%next
    enddo

  end subroutine med_fldList_GetFldNames

  !================================================================================

  subroutine med_fldList_Document_Mapping(logunit, med_coupling_active)

    ! input/output variables
    integer, intent(in)  :: logunit
    logical, intent(in)  :: med_coupling_active(:,:)

    ! local variables
    integer           :: nsrc,ndst
    integer           :: mapindex
    character(len=CS) :: mapnorm
    character(len=CL) :: mapfile
    character(len=CS) :: fldname
    character(len=CL) :: cvalue
    type(med_fldList_entry_type), pointer :: newfld
    character(len=*),parameter :: subname = '(med_fldList_Document_Mapping)'
    !-----------------------------------------------------------

    !---------------------------------------
    ! Document mapping (also add albedo and aoflux)
    !---------------------------------------

    ! Loop over src components
    do nsrc = 1,size(fldListFr)
       ! Loop over all possible destination components for each src component
       do ndst = 1,size(fldListTo)
          if (nsrc /= ndst .and. med_coupling_active(nsrc,ndst)) then
             ! Write all the mappings for fields from the src to the destination component
             write(logunit,*)' '
             newfld => fldListFr(nsrc)%fields
             do while(associated(newfld))
                mapindex = newfld%mapindex(ndst)
                if ( mapindex /= mapunset) then
                   call med_fld_GetFldInfo(newfld, compsrc=ndst, stdname=fldname, mapnorm=mapnorm, mapfile=mapfile, rc=rc)
                   if (chkerr(rc,__LINE__,u_FILE_u)) return

                   if (trim(mapnorm) == 'unset') then
                      cvalue = ' mapping '//trim(compname(nsrc))//'->'//trim(compname(ndst)) //' '//trim(fldname) // &
                           ' via '// trim(mapnames(mapindex))
                   else
                      cvalue = ' mapping '//trim(compname(nsrc))//'->'//trim(compname(ndst)) //' '//trim(fldname) // &
                           ' via '// trim(mapnames(mapindex)) // ' with '// trim(mapnorm) // ' normalization'
                   end if
                   write(logunit,100) trim(cvalue)
                   if (trim(mapfile) /= 'unset' .and. trim(mapfile) /= 'idmap') then
                      cvalue = ' and the mapping file '// trim(mapfile)
                      write(logunit,101) trim(cvalue)
                   end if
                end if
                newfld => newfld%next
             end do

          end if
       end do
    end do

    ! ocn-> atm mappings for atm/ocn fluxes computed in mediator on the ocn grid
    nsrc = compocn
    ndst = compatm
    if (med_coupling_active(nsrc,ndst) .and. allocated(fldlist_aoflux%fields%mapindex)) then
       newfld => fldlist_aoflux%fields
       do while(associated(newfld))
          call med_fld_GetFldInfo(newfld, compsrc=ndst, mapindex=mapindex, rc=rc)
          if ( mapindex /= mapunset) then
             call med_fld_GetFldInfo(newfld, stdname=fldname, compsrc=ndst, mapnorm=mapnorm, mapfile=mapfile, rc=rc)
             if (chkerr(rc,__LINE__,u_FILE_u)) return

             if (trim(mapnorm) == 'unset') then
                cvalue = ' mapping '//trim(compname(nsrc))//'->'//trim(compname(ndst)) //' '//trim(fldname) // &
                     ' via '// trim(mapnames(mapindex))
             else
                cvalue = ' mapping '//trim(compname(nsrc))//'->'//trim(compname(ndst)) //' '//trim(fldname) // &
                     ' via '// trim(mapnames(mapindex)) // ' with '// trim(mapnorm) // ' normalization'
             end if
             write(logunit,100) trim(cvalue)
             if (trim(mapfile) /= 'unset' .and. trim(mapfile) /= 'idmap') then
                cvalue = ' and the mapping file '// trim(mapfile)
                write(logunit,101) trim(cvalue)
             end if
          end if
          newfld => newfld%next
       end do
    end if

100 format(a)
101 format(3x,a)

  end subroutine med_fldList_Document_Mapping

  !================================================================================

  subroutine med_fldList_Document_Merging(logunit, med_coupling_active)

    !---------------------------------------
    ! Document merging to target destination fields
    !---------------------------------------

    ! input/output variables
    integer, intent(in)  :: logunit
    logical, intent(in)  :: med_coupling_active(:,:)

    ! local variables
    integer           :: nsrc,ndst
    character(len=CS) :: dst_comp
    character(len=CS) :: dst_field
    character(len=CS) :: src_comp
    character(len=CS) :: merge_type
    character(len=CS) :: merge_field
    character(len=CS) :: merge_frac
    character(len=CS) :: prefix
    character(len=CS) :: string
    character(len=CL) :: mrgstr
    type(med_fldList_entry_type), pointer :: newfld
    character(len=*),parameter :: subname = '(med_fldList_Document_Merging)'
    !-----------------------------------------------------------

    write(logunit,*)

    ! Loop over destination components
    do ndst = 1,size(fldListTo)
       dst_comp = trim(compname(ndst))
       prefix = '(merge_to_'//trim(dst_comp)//')'

       ! Loop over all flds in the destination component and determine merging data
       newfld => fldListTo(ndst)%fields
       do while(associated(newfld))
          call med_fld_GetFldInfo(newfld, stdname=dst_field, rc=rc)
          if (chkerr(rc,__LINE__,u_FILE_u)) return

          ! Loop over all possible source components for destination component field
          mrgstr = ' '
          do nsrc = 1,size(fldListFr)

             if (nsrc /= ndst .and. med_coupling_active(nsrc,ndst)) then
                src_comp    = compname(nsrc)
                call med_fld_GetFldInfo(newfld, compsrc=nsrc, merge_fields=merge_field, merge_type=merge_type, merge_fracname=merge_frac, rc=rc)
                if (chkerr(rc,__LINE__,u_FILE_u)) return

                if (merge_type == 'merge' .or. merge_type == 'sum_with_weights') then
                   string = trim(merge_frac)//'*'//trim(merge_field)//'('//trim(src_comp)//')'
                   if (mrgstr == ' ') then
                      mrgstr = trim(prefix)//": "// trim(dst_field) //'('//trim(dst_comp)//')'//' = '//trim(string)
                   else
                      mrgstr = trim(mrgstr) //' + '//trim(string)
                   end if
                else if (merge_type == 'sum') then
                   string = trim(merge_field)//'('//trim(src_comp)//')'
                   if (mrgstr == ' ') then
                      mrgstr = trim(prefix)//": "//trim(dst_field) //'('//trim(dst_comp)//')'//' = '//trim(string)
                   else
                      mrgstr = trim(mrgstr) //' + '//trim(string)
                   end if
                else
                   if (merge_type == 'copy') then
                      mrgstr = trim(prefix)//": " // trim(dst_field) //'('//trim(dst_comp)//')'//' = '// &
                           trim(merge_field)//'('//trim(src_comp)//')'
                   else if (merge_type == 'copy_with_weights') then
                      mrgstr = trim(prefix)//": "// trim(dst_field) //'('//trim(dst_comp)//')'//' = '// &
                           trim(merge_frac)//'*'//trim(merge_field)//'('//trim(src_comp)//')'
                   end if
                end if
             end if
          end do ! end loop over nsrc
          if (mrgstr /= ' ') then
             write(logunit,'(a)') trim(mrgstr)
          end if
          newfld => newfld%next
       end do ! end loop over fields
    end do  ! end loop over ndst

  end subroutine med_fldList_Document_Merging

end module esmflds
