module sh_recurrence
  use string,          only: str
  use iso_fortran_env, only: stdout => output_unit
  implicit none

  public :: sh_rec
  public :: sh_int
  public :: rec_type

  private

  real,parameter :: PI = 3.141592653589793

  type :: rec_type
     character            :: omiga= ''
     integer              :: l    = 0
     integer              :: m    = 0
     integer              :: num  = 0
     integer, allocatable :: index(:,:)  !< (2,num)
     real, allocatable    :: coef(:)     !< num
   contains
     procedure :: add_term => rec_add_term
     procedure :: clean   => rec_clean
     procedure :: print   => rec_print
  end type rec_type

contains
  !======================================================================
  subroutine rec_add_term(this,l,m,c)
    !< add term in recurrence relation
    class(rec_type), intent(inout) :: this
    integer, intent(in) :: l            !< index l
    integer, intent(in) :: m            !< index m
    real, intent(in)    :: c            !< recurrence coefficient

    integer, allocatable :: index_(:,:)
    real,    allocatable :: coef_(:)

    if (this%num > 0) then
       allocate(index_(1:2,1:this%num+1))
       allocate(coef_(1:this%num+1))

       index_(1:2,1:this%num) = this%index(1:2,1:this%num)
       coef_(1:this%num) = this%coef(1:this%num)

       this%num = this%num + 1
       call move_alloc(from=index_,to=this%index)
       call move_alloc(from=coef_,to=this%coef)
    else
       this%num = 1
       allocate(this%index(1:2, this%num))
       allocate(this%coef(this%num))
    end if

    this%index(1:2,this%num) = (/l,m/)
    this%coef(this%num) = c

  end subroutine rec_add_term


  !======================================================================
  subroutine rec_clean(this)
    !< clean recur
    class(rec_type), intent(inout) :: this

    if(allocated(this%index))  deallocate(this%index)
    if(allocated(this%coef))   deallocate(this%coef)

    this%l = 0
    this%m = 0
    this%num   = 0
    this%omiga = ''

  end subroutine rec_clean

  !======================================================================
  subroutine rec_print(this,unit)
    !< print recurrence relation
    class(rec_type), intent(in)    :: this
    integer,intent(in),optional :: unit
    integer :: lunit,i

    lunit = stdout
    if(present(unit)) lunit = unit

    write(lunit,'(a)') 'omiga('//this%omiga//')*'//'Y-'//str(this%l)//'-'//str(this%m)//' = &'
    do i = 1,this%num
       write(lunit,'(a)') str(this%coef(i))//'*Y-'//str(this%index(1,i))//'-'//str(this%index(2,i))
       if(i/= this%num) then
          write(lunit,'(a)') '+'
       end if

    end do


  end subroutine rec_print


  !======================================================================
  function sh_rec(omiga,n,m) result(this)
    !< add recurrence relation
    character, intent(in) :: omiga !< omiga = '1','x','y'
    integer, intent(in)   :: n     !< degree of spherical harmonics (n>=0)
    integer, intent(in)   :: m     !< order of spherical harmonics (-n<= m <= n)

    type(rec_type) :: this
    real        :: temp

    if( omiga /= '1' .and. omiga /= 'x'.and. omiga /= 'y'.and.omiga /= 'z') then
       error stop "error: The range of omiga: ['1','x','y','z']"
    end if

    if(n<0) then
       error stop 'error: The range of l: [l>=0]'
    end if

    if(abs(m)>n) then
       error stop 'error: The range of m: [-l<=m<=l]'
    end if

    this%omiga = omiga
    this%l     = n
    this%m     = m

    temp = 0.5/(2*n+1)
    if(omiga=='1')then
       call this%add_term(n,m,1.0)
       return
    elseif(omiga== 'x')then
       if(m>0)then
          call this%add_term(n-1,m+1,temp)
          call this%add_term(n+1,m+1,-temp)
          call this%add_term(n+1,m-1,temp*(n-m+1)*(n-m+2))
          call this%add_term(n-1,m-1,-temp*(n+m)*(n+m-1))
          return
       elseif(m==0)then
          call this%add_term(n-1,m+1,2*temp)
          call this%add_term(n+1,m+1,-2*temp)
          return
       elseif(m==-1)then
          call this%add_term(n-1,m-1,temp)
          call this%add_term(n+1,m-1,-temp)
          return
       elseif(m<-1)then
          call this%add_term(n-1,m-1,temp)
          call this%add_term(n+1,m-1,-temp)
          call this%add_term(n+1,m+1,temp*(n+m+1)*(n+m+2))
          call this%add_term(n-1,m+1,-temp*(n-m)*(n-m-1))
          return
       endif
    elseif(omiga== 'y')then
       if(m>1)then
          call this%add_term(n-1,-m-1,temp)
          call this%add_term(n+1,-m-1,-temp)
          call this%add_term(n-1,-m+1,-temp*(n-m+1)*(n-m+2))
          call this%add_term(n-1,-m+1,temp*(n+m)*(n+m-1))
          return
       elseif(m==1)then
          call this%add_term(n-1,-m-1,temp)
          call this%add_term(n+1,-m-1,-temp)
          return
       elseif(m==0)then
          call this%add_term(n-1,-1,2*temp)
          call this%add_term(n+1,-1,-2*temp)
          return
       elseif(m<0)then
          call this%add_term(n-1,-m+1,-temp)
          call this%add_term(n+1,-m+1,temp)
          call this%add_term(n+1,-m-1,temp*(n+m+1)*(n+m+2))
          call this%add_term(n-1,-m-1,-temp*(n-m)*(n-m-1))
          return
       endif
    elseif(omiga=='z')then
       call this%add_term(n+1,m,2*temp*(n-abs(m)+1))
       call this%add_term(n-1,m,2*temp*(n+abs(m)))
       return
    endif

  end function sh_rec

  !======================================================================
  function sh_int(omiga1,n1,m1,omiga2,n2,m2) result(res)
    !!< the integral of omiga1*sh1*omiga2*sh2
    character, intent(in) :: omiga1,omiga2 !< omiga = '1','x','y'
    integer, intent(in)   :: n1,n2         !< degree of spherical harmonics (n>=0)
    integer, intent(in)   :: m1,m2         !< order of spherical harmonics (-n<= m <= n)

    real        :: res
    type(rec_type) :: rec1, rec2
    integer     :: i,j

    rec1 = sh_rec(omiga1,n1,m1)
    rec2 = sh_rec(omiga2,n2,m2)

    res = 0.0
    do i = 1,rec1%num
       do j = 1,rec2%num
          if(rec1%index(1,i) == rec2%index(1,j).and.&
               rec1%index(2,i) == rec2%index(2,j)) then
             res = res + rec1%coef(i)*rec2%coef(j)/&
                  (sh_coef(rec1%index(1,i),rec1%index(2,i))*&
                  sh_coef(rec1%index(1,i),rec1%index(2,i)))
          end if
       end do
    end do

    res = res*sh_coef(n1,m1)*sh_coef(n2,m2)


  end function sh_int

  !======================================================================
  real function sh_coef(n,m)
    !< normalize coeffient
    integer,intent(in) :: n,m

    integer :: i
    real :: fact

    fact = 1.0d0

    if(m == 0) then
       sh_coef = sqrt((2.0*n+1.0)/(4.0*PI))
    else
       do i = n-abs(m)+1,n+abs(m)
          fact = fact*i
       end do
       sh_coef = sqrt((2.0*n+1.0)/(2.0*PI*fact))
    endif

  end function sh_coef

end module sh_recurrence
