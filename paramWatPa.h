!  paramWatPa.h - PORW31
!
      logical  linux7
      integer(C_INT) num_proc,kstart,np0,nq0,nr0,npqr0,    &
                     npqr30,npqr50,n00,npqio,mx,my,mz,mxy,mxyz, &
                     nbxs,ntmax,ip0,mintpol,nob3,iblk3
!
      character(len=29) praefixs,praefixi,praefixc,praefixe
      character(len=2)  suffix2,suffix1,suffix0
!     character   praefixs*29,praefixc*24,praefixe*24, &
!                 praefixi*24,suffix2*1
!
      parameter (num_proc=6,linux7=.true.)  ! 6 core, Linux
!     parameter (num_proc=8,linux7=.false.) 
!
      parameter (praefixs='/home/mtanaka/MPI_nano/PORW31')
      parameter (praefixi='/home/mtanaka/MPI_nano/porw31') ! read(12) - common by nfs
      parameter (praefixc='/home/mtanaka/MPI_nano/porw31') ! write(13)
      parameter (praefixe='/home/mtanaka/MPI_nano/porw31') ! WRITE_CONF
      parameter (suffix2='3a',suffix1='3a',suffix0='3',kstart=0)
!     parameter (suffix2='3b',suffix1='3a',suffix0='3',kstart=1) ! for restart
! 
      parameter (np0=50,nq0=300,nr0=95000)  !97000)  kl,km,kn=29x29x60, x 5 atoms
      parameter (npqr0=np0+nq0+nr0/5)       ! np+nq+nr/5
      parameter (npqr30=np0+nq0+3*nr0/5)    ! np+nq+3*nr/5
      parameter (npqr50=np0+nq0+nr0)        ! np+nq+nr
      parameter (npqio=np0+nq0+3*nr0/5)
!
      parameter (mx=81,my=81,mz=101) 
      parameter (mxy=mx*my,mxyz=mx*my*mz)   ! sendrecv
      parameter (n00=npqr30)
      parameter (nbxs=4000)
! 
      parameter (ntmax=7000) 
      parameter (ip0=3,mintpol=4*50048)     ! index number ip0
      parameter (nob3=13,iblk3=1)           ! nob3=13
