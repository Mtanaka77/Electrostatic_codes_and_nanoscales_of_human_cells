!  param_inv13.h
!
      integer(C_INT) &
                  np0,nq0,nq10,nq20,nr0,npq0,npqr0,MAXPART, &
                  kstart,maxmesh,mesh,IP,nhist,             &
                  mx,my,mz,mx1,my1,mz1,npio
!     character   praefixs*31,praefixc*31,praefixe*31, &
      character   praefixs*31,praefixc*24,praefixe*24, &
                  suffix2*2,suffix1*2,suffix0*1
!
      parameter  (np0=6,nq10=60+100,nq20=300,nr0=8000)
!
      parameter  (nq0=nq10+nq20,nhist=10000)
      parameter  (npqr0=np0+nq0+nr0,MAXPART=npqr0)
      parameter  (npq0=np0+nq0,npio=np0+nq0)
!
      parameter  (maxmesh=32,mesh=32,IP=3)
      parameter  (mx=mesh,my=mesh,mz=mesh)
      parameter  (mx1=mesh+1,my1=mesh+1,mz1=mesh+1)
!
      parameter  (praefixs='/home/tanakam/MPI_chinv3/CIMV13', & 
                  praefixc='/data/sht/tanakam/cimv13',        &
                  praefixe='/data/sht/tanakam/CIMV13')
!     parameter  (praefixs='/home/mtanaka/MPI_chinv3/CIMV23', & 
!                 praefixc='/home/mtanaka/MPI_chinv3/cimv23', &
!                 praefixe='/home/mtanaka/MPI_chinv3/CIMV23')
!
      parameter  (kstart=0,suffix2='3a',suffix1='3a',suffix0='3')
!
      integer(C_INT) P_max,MINTPOL,Brillouin
      parameter  (P_max=3,MINTPOL=4*50048,Brillouin=1)
!
      real(C_DOUBLE) A1,A2,A3,A4,A5,PP
      parameter ( A1 = 0.254829592, A2 = -0.284496736 )
      parameter ( A3 = 1.421413741, A4 = -1.453152027 )
      parameter ( A5 = 1.061405429, PP  =  0.3275911  )

