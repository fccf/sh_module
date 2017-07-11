program test_shdof
  use sh_dof
  implicit none

  type(shdof) :: sd

  sd = shdof(pnmax=5,dim=2)

  call sd%print()

end program test_shdof
