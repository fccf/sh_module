module sh_dof
  use iso_fortran_env, only: error_unit,output_unit
  implicit none

  public :: shdof

  private

  !!>the degree of freedom of spherical harmonics
  type :: shdof
     !!> index type
     integer :: type  = 1
     !!> max order of spherical harmonic
     integer :: pnmax = 0
     !!> number of dimension
     integer :: dim   = 0
     !!> number of freedom 1d:pn+1, 2d:(pn+1)(pn+2)*0.5, 3d:(pn+1)*(pn+1),
     integer :: nm    = 0
     !!> max degree of freedom md = (pn+1)*(pn+1)
     integer :: md    = 0
     !!> the nmj table (md,3)
     integer, allocatable :: nmj(:,:)
     !!> dof (nm)
     integer, allocatable :: dof(:)
   contains
     procedure :: clear
     procedure :: print
     procedure :: index_j
     procedure :: index_n
     procedure :: index_nm

  end type shdof

  interface shdof
     module procedure init
  end interface shdof

  integer, parameter :: SHDOF_TYPE1 = 1
  integer, parameter :: SHDOF_TYPE2 = 2


contains
  !======================================================================
  function init(pnmax,dim,type) result(this)

    integer, intent(in) :: pnmax
    integer, intent(in), optional :: dim
    integer, intent(in), optional :: type

    type(shdof) :: this

    integer :: ldim, ltype

    if(present(dim)) then
       ldim = dim
    else
       ldim = 3
    end if

    if(present(type)) then
       ltype = type
    else
       ltype = SHDOF_TYPE1
    end if

    if (ltype==SHDOF_TYPE1 .and. ltype==SHDOF_TYPE2) then
       error stop 'error: The type of shdof is wrong! the type should be 1 or 2.'
    end if

    call this%clear()

    this%dim   = ldim
    this%pnmax = pnmax
    this%type  = ltype
    this%md = (pnmax+1)*(pnmax+1)

    if(this%dim == 1) then
       this%nm = pnmax+1
    elseif(this%dim == 2) then
       this%nm = (pnmax+1)*(pnmax+2)*0.5
    elseif(this%dim == 3)then
       this%nm = this%md
    else
       error stop 'error: dimension error!'
    endif

    allocate(this%nmj(this%md,3))
    allocate(this%dof(this%nm))
    this%nmj = 0
    this%dof = 0

    call cal_nmj(this)
    call set_dof(this)

  end function init
  !======================================================================
  subroutine clear(this)

    class(shdof), intent(in out) :: this

    this%type  = 0
    this%dim   = 0
    this%pnmax = 0
    this%nm    = 0

    if(allocated(this%nmj))   deallocate(this%nmj)
    if(allocated(this%dof))   deallocate(this%dof)

  end subroutine clear
  !======================================================================
  subroutine print(this,unit)

    class(shdof),intent(inout)   :: this
    integer,intent(in), optional :: unit

    integer :: i,j, lunit

    lunit = output_unit
    if(present(unit)) lunit = unit


    write(lunit,'(1x,a)') 'The DOF information of spherical harmonics'
    write(lunit,'(1x,a,i3)') 'The index type                    =',this%type
    write(lunit,'(1x,a)')    'nmjtable                          = ...'
    write(lunit,'(1x,3i15)') ((this%nmj(i,j),j=1,3),i=1,this%md)

    write(lunit,'(1x,a,i3)') 'Number of space dimension         = ',this%dim
    write(lunit,'(1x,a,i3)') 'The order of spherical harmonics  = ',this%pnmax
    write(lunit,'(1x,a)')    'dofarray                          = ...'
    write(lunit,'(1x,i10)')(this%dof(i),i=1,this%nm)

  end subroutine print
  !======================================================================
  subroutine cal_nmj(this)

    class(shdof), intent(inout) :: this
    integer :: n,m,j

    do n = 0,this%pnmax
       do m = -n,n

          j = this%index_j(n,m)
          this%nmj(j,1) = n
          this%nmj(j,2) = m
          this%nmj(j,3) = j

       end do
    end do

  end subroutine cal_nmj
  !======================================================================
  subroutine set_dof(this)

    class(shdof),intent(in out) :: this

    integer :: n,m,j,j1,j2,j3,nm(2)

    j1 = 0;j2=0;j3=0

    do j = 1,this%md

       nm = this%index_nm(j)
       n = nm(1)
       m = nm(2)
       if(this%dim == 1) then
          if(m == 0) then
             j1 = j1+1
             this%dof(j1) = j
          endif
       elseif(this%dim == 2)then
          if(mod(m+n,2)==0) then  !XY plane  2017/1/14
          !if( m ==0.or.(m>0.and.mod(m,2)==0).or.(m<0.and.mod(m,2)/=0)) then !YZplane
          !if(m>=0) then           !ZX plane
             j2 = j2+1
             this%dof(j2) = j
          endif
       elseif(this%dim == 3) then
          j3 = j3+1
          this%dof(j3) = j
       endif

    end do


  end subroutine set_dof
  !======================================================================
  function index_j(this,n,m) result(j)

    class(shdof), intent(in) :: this
    integer,intent(in) :: n,m
    integer :: j

    j = 0
    if(this%type == SHDOF_TYPE1) then
       j = n*n+n+m+1
    elseif(this%type == SHDOF_TYPE2) then
       if( m>=0 ) then
          j = n*n + 2*m + 1
       else
          j = n*n -2*m
       end if
    end if

  end function index_j
  !======================================================================
  function index_n(this,j) result(n)

    class(shdof), intent(in) :: this
    integer,intent(in) :: j
    integer :: n

    n = this%nmj(j,1)

  end function index_n
  !======================================================================
  function index_nm(this,j) result(nm)

    class(shdof), intent(in) :: this
    integer,intent(in) :: j
    integer :: nm(2)

    nm(1) = this%nmj(j,1)
    nm(2) = this%nmj(j,2)

  end function index_nm
  !

end module sh_dof
