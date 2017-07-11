module sh_module

  use sh_value
  use sh_recurrence
  use sh_dof
  use petsc_matrix

  implicit none

  public :: sh_init
  public :: sh_clear
  public :: GAX,GAY,GAZ,GAN,GSHDOF

  private

  type(petsc_mat)  :: GAX,GAY,GAZ,GAN
  type(shdof)      :: GSHDOF

  real,parameter :: EPS_ZERO = 1.0e-12

contains
  !======================================================================
  subroutine SH_init(pn,nd,tp)

    integer, intent(in) :: pn,nd,tp

    GSHDOF = shdof(pn,nd,tp)

    call cal_Am()

  end subroutine SH_init
  !======================================================================
  subroutine SH_clear()

    call GSHDOF %clear()
    call GAX%clear()
    call GAY%clear()
    call GAZ%clear()

  end subroutine SH_clear
  !======================================================================
  subroutine cal_Am()

    integer :: i,j,l1,m1,l2,m2,d1,d2,nm1(2),nm2(2)
    real    :: v

    call GAX%creat(GSHDOF%nm,GSHDOF%nm)
    call GAY%creat(GSHDOF%nm,GSHDOF%nm)
    call GAZ%creat(GSHDOF%nm,GSHDOF%nm)

    v = 0.0
    do i = 1,GSHDOF%nm
      d1 = GSHDOF%dof(i)
      nm1 = GSHDOF%index_nm(d1)
      l1 = nm1(1);m1 = nm1(2)

      do j = 1,GSHDOF%nm
        d2 = GSHDOF%dof(j)
        nm2 = GSHDOF%index_nm(d2)
        l2 = nm2(1);m2 = nm2(2)

        v  =sh_int('1',l1,m1,'x',l2,m2)
        if(abs(v)>EPS_ZERO) then
          call GAX%set(i-1,j-1,v)
        endif

        v  =sh_int('1',l1,m1,'y',l2,m2)
        if(abs(v)>EPS_ZERO) then
          call GAY%set(i-1,j-1,v)
        endif

        v  =sh_int('1',l1,m1,'z',l2,m2)
        if(abs(v)>EPS_ZERO) then
          call GAZ%set(i-1,j-1,v)
        endif

      end do
    end do

    call GAX%assemble()
    call GAY%assemble()
    call GAZ%assemble()

  end subroutine cal_Am

end module sh_module
