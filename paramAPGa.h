!  paramAPGa.h - PORW21
!
      logical  linux7
      integer(C_INT) num_proc,kstart,np0,nq0,nr0,npq0,npqr0, &
                     npio,mx,my,mz,mxy,mxyz,nbxs,ntmax,      &
                     ip0,mintpol,nob3,iblk3
!
      character(len=29) praefixs,praefixi,praefixc,praefixe
      character(len=2)  suffix2,suffix1,suffix0
!     character   praefixs*29,praefixc*24,praefixe*24, &
!                 praefixi*24,suffix2*1
!
      parameter (num_proc=6,linux7=.true.)  ! 6 core, Linux
!     parameter (num_proc=8,linux7=.false.) 
!
      parameter (praefixs='/home/mtanaka/MPI_nano/PORW21')
      parameter (praefixi='/home/mtanaka/MPI_nano/porw21') ! read(12) - common by nfs
      parameter (praefixc='/home/mtanaka/MPI_nano/porw21') ! write(13)
      parameter (praefixe='/home/mtanaka/MPI_nano/porw21') ! WRITE_CONF
      parameter (suffix2='3a',suffix1='3a',suffix0='3',kstart=0)
!     parameter (suffix2='3b',suffix1='3a',suffix0='3',kstart=1) ! for restart
!
      parameter (np0=50,nq0=300,nr0=15000)  ! kl,km,kn=27x27x57
      parameter (npq0=np0+nq0,npqr0=np0+nq0+nr0)
      parameter (npio=np0+nq0)
!
      parameter (mx=81,my=81,mz=101) 
      parameter (mxy=mx*my,mxyz=mx*my*mz)   ! sendrecv
      parameter (nbxs=3700) ! 3200 ok
!
      parameter (ntmax=7000) 
      parameter (ip0=3,mintpol=4*50048)     ! index number ip0
      parameter (nob3=13,iblk3=1)           ! nob3=13
