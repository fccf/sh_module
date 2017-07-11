program test_shmodule
  use sh_module
  use fu_module
  use petsc_interface
  implicit none

  type(petsc)        :: PS
  type(petsc_mat)    :: DM, MOM, RM

  integer :: i,j
  integer :: nd  = 3
  integer :: pn = 1
  integer :: tp = 1
  integer :: nn = 20

  real,pointer :: rmx(:,:)

  call PS%init()
  ! call get_option('-pn',pn)
  ! call get_option('-nd',nd)
  ! call get_option('-tp',tp)

  call SH_init(pn,nd,tp)

  call gAx%print('Ax')
  call gAx%toDense(DM)
  call DM%print('Axm')
  call DM%clear()

  call gAy%print('Ay')
  call gAy%toDense(DM)
  call DM%print('Aym')
  call DM%clear()

  call gAz%print('Az')
  call gAz%toDense(DM)
  call DM%print('Azm')
  call DM%clear()


  call RM%creat(nn,gShdof%nm,'DENSE')
  call MOM%creat(nn,gShdof%nm,'DENSE')
  do i = 1,MOM%nrow
    do j = 1,MOM%ncol
      call MOM%set(i-1,j-1,1.0)
    enddo
  enddo
  call MOM%assemble()

  call MOM%print('MMT')

  call MOM%mm(gAx,RM)
  call RM%print('RMMX')

  call RM%get(rmx)

  write(*,*) 'Size RMX:',size(rmx,1),size(rmx,2)
  do i = 1,size(rmx,1)
    write(*,*) (rmx(i,j),j=1,size(rmx,2))
  enddo
  rmx(:,2) = 22222222.
  call RM%restore(rmx)
  call RM%print('RMMX1')

  call MOM%mm(gAy,RM)
  call RM%print('RMMY')
  call MOM%mm(gAz,RM)
  call RM%print('RMMZ')

  call PS%clear()



end program test_shmodule
