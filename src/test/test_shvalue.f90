program test_shvalue

  use sh_value
  implicit none

  call testPl()
  call testPlm()

contains
  !======================================================================
  subroutine testPl()

    integer,parameter :: NN   = 100
    integer,parameter :: lmax = 5

    real,parameter    :: PI = 3.1415926535897
    real :: pl(lmax+1),x
    real :: theta(NN)
    integer :: i,l,li,pl_unit

    open(newunit = pl_unit,file='../data/pl.dat')

    do i = 1,NN
       theta(i) = i*PI/NN
       x = cos(theta(i))
       li = 0
       do l = 0,lmax
          li = li+1
          pl(li) = Plm(l,0,x)
       end do

       write(pl_unit,'(100es15.7)') x,pl

    end do

    close(pl_unit)

  end subroutine testPl
  !======================================================================
  subroutine testPlm()
    integer,parameter :: NN   = 100
    integer,parameter :: lmax = 5

    real,parameter    :: PI = 3.1415926535897
    real :: pnm((lmax+1)*(lmax+2)/2),x
    real :: theta(NN)
    integer :: i,l,m,lmi
    integer :: plm_unit

    open(newunit = plm_unit,file='../data/plm.dat')

    do i = 1,NN
       theta(i) = i*PI/NN
       x = cos(theta(i))
       lmi= 0
       do l = 0,lmax
          do m = 0,l
             lmi = lmi+1
             pnm(lmi) = Plm(l,m,x)*alm(l,m)*sqrt(4.0*PI/(2.0*l+1.0))

          end do

       end do

       write(plm_unit,'(100es15.7)') x,pnm

    end do

    close(plm_unit)

  end subroutine testPlm



end program test_shvalue
