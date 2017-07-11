module sh_value
  implicit none

  public :: ylm,plm,alm

  private

  real,parameter :: PI = 3.141592653589793


contains
  !======================================================================
  real function ylm(l,m,theta,phi)
    !!< Calculate the value of real spherical harmonics

    integer,intent(in) :: l      !< Degree of Ylm
    integer,intent(in) :: m      !Order of Ylm
    real, intent(in)   :: theta  !<the latitudinal coordinates
    real, intent(in)   :: phi    !the longitudinal coordinates

    real :: coef

    coef = alm(l,m)

    if( m>=0 ) then
       Ylm=coef*Plm(l,m,cos(theta))*cos(m*phi)
       return
    else
       Ylm=coef*Plm(l,-m,cos(theta))*sin(-m*phi)
    end if

  end function ylm
  !======================================================================
  real function alm(l,m)

    integer,intent(in) :: l,m

    integer :: i
    real :: fact

    fact = 1.0

    if(m == 0) then
       alm = sqrt((2.0*l+1.0)/(4.0*PI))
    else
       do i = l-abs(m)+1,l+abs(m)
          fact = fact*i
       end do
       alm = sqrt((2.0*l+1.0)/(2.0*PI*fact))
    endif

  end function alm

  !======================================================================
  real function plm (l,m,x)
    !!< Calculate the value of Associated Legendre Polynomial by recursion relation
    integer,intent(in) :: l, m
    real, intent(in) :: x

    real :: Pmm, Pmp1m, fact, x2
    integer :: i

    Plm = 0.0
    Pmm = 1.0
    if(m>=0) then
       fact = 1.0
       x2 = sqrt(1.0-x*x)
       do i=1,m
          Pmm = Pmm*(-fact)*x2
          fact=fact + 2.0
       end do
    endif

    Pmp1m=x*(2*m+1)*Pmm

    if(l==m)then
       Plm=Pmm
    elseif(l==m+1)then
       Plm=Pmp1m
    elseif( l>=m+2 ) then
       do i=m+2,l
          Plm=((2.0*i-1.0)*x*Pmp1m-(i+m-1.0)*Pmm)/(i-m)
          Pmm=Pmp1m
          Pmp1m=plm
       end do
    end if

  end function plm


end module sh_value
