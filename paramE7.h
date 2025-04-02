! paramE7.h
      integer(C_INT) num_proc,np0,nq0,npq0,MAXPART
      real(C_DOUBLE) R_sp0,D_sp0
      character  praefixc*6,cname*1
!
      parameter  (num_proc=6)  ! Linux
!     parameter  (num_para=8)  ! NIFS
      parameter  (praefixc='expl07',cname='a')
!
      parameter  (np0=4200,nq0=4200)
      parameter  (R_sp0=1.0d-6,D_sp0=1.0d21)
      parameter  (npq0=np0+nq0,MAXPART=np0+nq0)

