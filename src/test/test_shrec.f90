program test_shrec
  use sh_recurrence
  use sh_dof
  use fu_module
  implicit none

  type(rec_type) :: REC
  type(shdof) :: SDof
  integer :: lmax = 1,tp = 1
  integer :: d
  real,allocatable :: a1(:,:),a1x(:,:),a1y(:,:),a1z(:,:)
  integer :: i,j,l1,m1,l2,m2,nm1(2),nm2(2)

  REC = sh_rec('x',4,2)
  call REC%print()

  call get_Option('-pn',lmax)
  call get_Option('-tp',tp)

  SDof = shdof(lmax,3,tp)

  d = SDof%nm
  allocate(a1(d,d),a1x(d,d),a1y(d,d),a1z(d,d))

  do i = 1,SDof%nm
    nm1 = SDof%index_nm(i)
    l1 = nm1(1);m1 = nm1(2)
    do j = 1,SDof%nm
      nm2 = SDof%index_nm(j)
      l2 = nm2(1);m2 = nm2(2)
      a1(j,i) = sh_int('1',l1,m1,'1',l2,m2)
      a1x(j,i) = sh_int('1',l1,m1,'x',l2,m2)
      a1y(j,i) = sh_int('1',l1,m1,'y',l2,m2)
      a1z(j,i) = sh_int('1',l1,m1,'z',l2,m2)
    end do
  end do

  call save('a1',a1)
  call save('a1x',a1x)
  call save('a1y',a1y)
  call save('a1z',a1z)


end program test_shrec
