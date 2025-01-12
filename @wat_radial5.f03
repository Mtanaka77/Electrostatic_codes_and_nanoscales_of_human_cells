!**************************************************************
!*  Snapshots of radial distribution functions.               *
!*                                 1/24/2020 Motohiko Tanaka  *
!*  gfortran @wat_radial5.f03                                 *
!*  pgf95 @wat_radial5.f03                                    *
!***************************** 10/14/2005 **** 2/06/2018 ******
! 
!  << Remark >> write FT13 in real*4
!
      program wat_redial
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer(C_int) np0,nq10,nq20,nr0,npq0,MAXPART, &
                     maxmesh,mesh,IP,mx,my,mz,mx1,my1,mz1
!
!     parameter  (np0=6,nq10=60+100,nq20=300) 
!     parameter  (np0=3,nq10=30+100,nq20=300) 
      parameter  (np0=1,nq10=10+100,nq20=300) 
!
      parameter  (npq0=np0+nq10+nq20,MAXPART=npq0)
      parameter  (maxmesh=32,mesh=32,IP=3)
      parameter  (mx=mesh,my=mesh,mz=mesh)
      parameter  (mx1=mesh+1,my1=mesh+1,mz1=mesh+1)
!
      character     fig_label*32
!
      real(C_float),dimension(npq0) :: &
                    x4,y4,z4,vx4,vy4,vz4,ch4,am4,ag4
      real(C_float) t_unit,a_unit,w_unit,e_unit,pi2,tequil,econv,gamma, &
                    Qcore4,Rmac4,Wmac4,Zcp4,Zcn4,qfrac4,alpha4,edc4,    &
                    tau_wave4,xmax4,ymax4,zmax4,t0(100),tmax,tmin,      &
                    hh,exc,a1,a2,t,xleng,tau_wave,bjerrum4,ep4,dd
!
      integer(C_int)  np,nq,npio,nq1,nq2,nseg,nps,ic,icmax,nplot,npq, &
                    avsplo(npq0),nsg(30),ipl,iffig,           &
                    i,j,k,n1,n2,n3,nframe,iplot,is,it
!
      character     praefixs*6,praefixo*6,postfix(20)*2,     &
                    midfix*6,suffix*3,isymb*4,               &
                    LABEL,cdate*10,ctime*8
      logical       first
!
      common/parm1/ xmax4,ymax4,zmax4
      common/parmd/ cdate,ctime
      COMMON/HEADR2/ t,xleng,Edc4,tau_wave
!
!     namelist/inp1/ tmin,tmax
!
      call date_and_time_7 (cdate,ctime)
      write(6,*) 'Today = ',cdate,'  time = ',ctime
!
!--------------------------------------------------
!*  pgf95 @wat_radial5.f03   
!*  gfortran @wat_radial5.f03      
!
      tmin=      0.
      tmax= 110001. !<--
      tequil=    0. !<--
!
      midfix='cimv15'  !<--
      icmax= 1         !<-- 3a
!
      postfix(1)= '3a' !<--
!     postfix(2)= '3b' 
!
      iffig= 0         ! =1 for no labels
      nplot= 70
!
      first  = .true.
      suffix='.13'
!
      do 10 i= 1,nplot
      t0(i)= 500.*i ! 50ps(?) interval, 1/2 cycles
   10 continue             !       ****
!
      tmax= amin1(tmax,t0(nplot)+1.)
!
!     write(6,*) '&inp1 tmin,tmax ...'
!--------------------------------------------------
!
      fig_label= 'CHG_INV '
      praefixs= midfix ! cimv21
      praefixo= midfix ! 
!
      it= 0
      is= 0
      ipl= 0
!
      nframe= 4
      OPEN (unit=77,file=midfix//'_78.ps')
      write(6,*) ' Write: ',midfix//'_78.ps'
!
      call gopen (nframe)
!
!
      ic= 1
      write(6,*) 'Read: ',midfix//'.13'//postfix(1)
      OPEN (unit=13,file=midfix//'.13'//postfix(1),form='unformatted') 
!
    1 continue
      read(13) np,nq,npio,nseg,qfrac4,Qcore4,Rmac4,Wmac4,     &
               Zcp4,Zcn4,Bjerrum4,alpha4,edc4,xmax4,ymax4,zmax4
      read(13) ch4,am4,ag4,ep4
!
      nq1= nq10
      nq2= nq20
      npq= np +nq
!
      n1= 0
      n2= 0
!
      write(6,50) np,nq1,nq2,npq
   50 format('np,nq1,nq2,npq=',4i7)
!
!  gfortran @wat_radial5.f03
!     Maximum wave number 20= dd*50 -> dd=20./50
      dd=  20./50
!
      if(ic.eq.1) then
        write(6,*) 'Rmac4, dd=',Rmac4,dd
      end if
!
      do i= np+1,npq
      if(ch4(i).ge.0.d0) then
        n1= n1 +1
        a1= ag4(i)
      end if
      if(ch4(i).lt.0.) then
        n2= n2 +1
        a2= ag4(i)
      end if
      end do
!  -------------------------
!
      nseg= 0
      if(nseg.ne.0) then
        do 30 k= 1,nseg+1
        nsg(k)= 1 + nps*(k-1)
   30   continue
!
        if(mod(n1,nps).ne.0) nsg(nseg+1)= n1 +1
        write(6,*) ' nps, nseg, last i=',nps,nseg,nsg(nseg+1)
      end if
!
      if(ic.eq.1) then
        write(6,*) ' Number of: np, N(+), N(-)=',np,nq1,nq2
        write(6,*) ' ion radius: Rmac4, a(+), a(-)=',Rmac4,a1,a2
        write(6,*) ' valence: Qcore, Zcp, Zcn=',Qcore4,Zcp4,Zcn4
        write(6,*) ' Edc4=',edc4
      end if
!
!     t_unit= 0.01d-12          ! 0.01 ps
!     a_unit= 1.00d-08          ! 1 Ang
!     w_unit= 1.6605d-24
!     e_unit= 4.803d-10
      pi2   = 8*datan(1.0d0)
!
!     econv = t_unit**2*e_unit/(w_unit*a_unit*300) ! V/cm
!     edc4= econv * edc4
!
!*  gfortran @wat_radial5.f03
  100 read(13,end=3) t,x4,y4,z4
      it= it +1
!
      if(t.ge.tequil) then
        exc = Edc4 
      else
        exc= 0
      end if
!
!-------------------------------------
!* Avoid overlap (Restart failure)
!-------------------------------------
!
      if(t.lt.tmin) go to 100
      if(t.gt.tmax) go to 1000
!
      if(mod(it,100).ne.1) goto 100 
      is= is +1
!       write(6,*) 'is,t=',is,t
!
      iplot= 0
      do k= 1,nplot
      if(abs(t-t0(k)).lt.1.0) then
        iplot= 1
        ipl= ipl +1
!
        hh= 0.7
        call symbol (0.0,18.0,hh,midfix//'.13'//postfix(1),0.,11)
        write(6,*) 'it, t= ',it,t
        go to 70
      end if
      end do
      go to 100
!
! --------------------------------
!*  Radial distribution functions
! --------------------------------
   70 continue
      call rdistr (x4,y4,z4,ch4,Zcp4,Zcn4,Rmac4,dd,  &
                                     np,nq1,nq2,iffig)
      go to 100
!!
    3 ic= ic +1
      if(ic.gt.icmax) go to 1000
!
      OPEN (unit=13,file=midfix//'.13'//postfix(ic),form='unformatted') 
        write(6,*) 'ic, FT13=',ic,midfix//'.13'//postfix(ic)
      go to 1
!
 1000 call gclose
      close(13)
!
      write(6,*)
      write(6,*) ' gfortran @wat_radial5.f03'
      write(6,*) ' write to: ',midfix//'_78.ps'
      write(6,*) '  then pdf to: pspdf ',midfix//'_78(.pdf)'
!
      stop
      end program wat_redial
!
!
!------------------------------------------------------------------------
      subroutine date_and_time_7 (date_now,time_now)
!------------------------------------------------------------------------
      implicit none
!
      integer, dimension(8) :: ipresent_time
      character(len=8)  :: time_now
      character(len=10) :: date_now

      call date_and_time(values=ipresent_time)

      write(time_now,'(i2,":",i2,":",i2)') ipresent_time(5:7)
      write(date_now,'(i4,"/",i2,"/",i2)') &
               ipresent_time(1),ipresent_time(2),ipresent_time(3)
!
      return
      end subroutine date_and_time_7
!
!
!-------------------------------------------------------------------
      subroutine rdistr (x,y,z,ch,Zcp,Zcn,Rmac4,dd,np,nq1,nq2,iffig)
!-------------------------------------------------------------------
      use, intrinsic :: iso_c_binding 
      implicit none
!
      integer*4  np0,nq10,nq20,npq0,kmax
!
!     parameter  (np0=6,nq10=60+100,nq20=300)
!     parameter  (np0=3,nq10=30+100,nq20=300)
      parameter  (np0=1,nq10=10+100,nq20=300)
      parameter  (kmax=51,npq0=np0+nq10+nq20)
!
      real(C_float),dimension(npq0) :: x,y,z,ch
      real(C_float),dimension(kmax) :: hh,accum,hmpn,hmp,hmn
      real(C_float) Zcp,Zcn,f1max,f2max,f3max,f4max,      &
                    f1min,f2min,f3min,f4min,chmac,        &
                    pi,dd,Rmac4,dx,dy,dz,rr,f7max,f7min
!
      integer(C_int) np,nq,nq1,nq2,npq,iffig,iffig2,imx,  &
                     i,j,k,ir,nxmax,ILN,ILG,nxtick,nytick
!
      real(C_float) xmax,ymax,zmax,time,xleng,Edc4,tau_wave
      character  cdate*10,ctime*8
      common/parm1/ xmax,ymax,zmax
      common/parmd/ cdate,ctime
      COMMON/HEADR2/ time,xleng,Edc4,tau_wave
!
      np = np0
      nq = nq1 +nq2
      npq= np +nq1 +nq2
!
      pi= 4.*atan(1.0)
!
      do k= 1,kmax
      hmp(k) = 0
      hmn(k) = 0
      hmpn(k)= 0
!
      hh(k)= dd*(k-1)
      end do
!
!* gfortran @wat_radial5.f03
      chmac= 0
      do i= 1,np
      chmac= chmac +ch(i)
      end do
      chmac= chmac/np     ! per macroion
!
      imx= Rmac4/dd
      do i= 1,imx
      hmpn(i)= chmac/imx  ! macroion is scattered at i=1,imx
      end do
!
      accum(1)= chmac     ! macroion stands at r= 0
!
!* Disribution of counterions and coions
      do i= 1,np
      do j= np+1,npq
!
      dx= x(i) -x(j)
      dy= y(i) -y(j)
      dz= z(i) -z(j)
!
      dx = dx - nint(dx/xmax)*xmax
      dy = dy - nint(dy/ymax)*ymax
      dz = dz - nint(dz/zmax)*zmax
      rr= sqrt(dx**2 +dy**2 +dz**2) 
!
!   max number 20= dd*50 -> dd=20./50
      ir= rr/dd +1.001
!
      if(ir.ge.1 .and. ir.le.kmax) then
        if(ch(j).gt.0.d0) then
          hmp(ir) = hmp(ir) +ch(j)/np  ! per macroion
!
        else if(ch(j).lt.0.d0) then
          hmn(ir) = hmn(ir) +abs(ch(j))/np
        end if
!                                  macroion at the origin is not counted !
        hmpn(ir)= hmpn(ir) +ch(j)/np  !! positive and negative
      end if
      end do
      end do
!
!
      do i= 2,kmax
      accum(i)= accum(i-1) +hmpn(i)
      end do
!
      CALL SYMBOL ( 5.0,18.0,0.7,'@radial5.f03',0.,12)
      CALL SYMBOL ( 9.0,18.0,0.7,'date: ',0.,6)
      CALL SYMBOL (10.6,18.0,0.7,cdate,0.,10)
      CALL SYMBOL (14.5,18.0,0.7,'time: ',0.,6)
      CALL SYMBOL (16.1,18.0,0.7,ctime,0.,8)
!
      CALL SYMBOL (11.0,17.2,0.7,'Edc=',0.,4)
      CALL VALUES (12.5,17.2,0.7,edc4,0.,101)
!     CALL SYMBOL (17.0,17.2,0.7,'tau=',0.,4)
!     CALL VALUES (18.5,17.2,0.7,tau_wave,0.,101)
!
      nxmax= kmax ! just cut
!
      call lplmax (hmp,  f1max,f1min,nxmax)
      call lplmax (hmn,  f2max,f2min,nxmax)
      f7max= amax1(f1max,f2max)
      f7min= amin1(f1min,f2min)
!
      call lplmax (hmpn, f3max,f3min,nxmax)
      call lplmax (accum,f4max,f4min,nxmax)
!       write(6,*) 'f3max,f3min=',f3max,f3min
!       write(6,*) 'f7max,f7min=',f7max,f7min
!       do i= 1,kmax
!       write(6,301) i,hmpn(i),accum(i),hmp(i),hmn(i)
! 301   format(i5,1p4e12.3)
!       end do
!
      ILN= -1
      ILG= 2
! 
      nxtick= 9
      nytick= 5
!       fmax7= 2.3e-2 ! Ookiku suru to Shukusou !!
!       +++++++++++++
!
      call hplot1 (2,2,nxmax,hh, hmpn,f3max,f3min,ILN,nxtick,nytick, &
                    'M-co/cou',8,'        ',8,'        ',8,iffig)
      call hplot1 (2,3,nxmax,hh,accum,f4max,f4min,ILN,nxtick,nytick, &
                    'Accum fn',8,'        ',8,'        ',8,iffig)
!
      call hplot1 (3,2,nxmax,hh,  hmp,f7max,0.,   ILN,nxtick,nytick, &
                    'counter ',8,'        ',8,'        ',8,iffig)
      call hplot1 (3,3,nxmax,hh,  hmn,f7max,0.,   ILN,nxtick,nytick, &
                    'coion   ',8,'        ',8,'        ',8,iffig)
      call chart
!
      return
      end subroutine rdistr
!
!
!------------------------------------------------------
      subroutine lplmax (f,fmax,fmin,is)
!------------------------------------------------------
      use, intrinsic :: iso_c_binding  ! <-
      implicit none 
!
      integer(C_int) i,is
      real(C_float),dimension(is) :: f
      real(C_float)  fmax,fmin
!
      fmax= -1.e10
      fmin=  1.e10
!
      do i= 1,is
      fmax= amax1(fmax,f(i))
      fmin= amin1(fmin,f(i))
      end do
!
      return
      end subroutine lplmax
!
!
!-------------------------------------------------
       subroutine circle (x,y,d,ic)
!-------------------------------------------------
!*  Open circle centered at (x,y) /or outer edge.
!
      write(77,*) " 3.0 setlinewidth"
!
      pi= 3.1415927
      nc= 13
      dth= 2.*pi/nc
      a= d/2.
!
      x0= x +a
      y0= y
      call plot (x0,y0,3)
!
      do 100 j= 1,nc
      th= dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      call plot (x1,y1,2)
  100 continue
!
      call plot (x1,y1,3)
      write(77,*) " 1.0 setlinewidth"
!
      if(ic.eq.1) return
!------------------------------------
!*  Filled circle centered at (x,y).
!------------------------------------
!
      write(77,*) " 3.0 setlinewidth"
!
      nc= 18
      dth= pi/(2*nc +1)
!
      do 300 j= -nc,nc
      th= 0.5*pi +dth*j
!
      x1= x +a*cos(th)
      y1= y +a*sin(th)
!
      x2= x1
      y2= 2.*y -y1
!
      call plot (x1,y1,3)
      call plot (x2,y2,2)
  300 continue
!
      call plot (x2,y2,3)
      write(77,*) " 3.0 setlinewidth"
      write(77,*) " 1.0 setlinewidth"
!
      return
      end
!
!
!-------------------------------------------------
        subroutine triang (x,y,d,ic)
!-------------------------------------------------
!*  Open triangle centered at (x,y) /or outer edge.
!
      write(77,*) " 3.0 setlinewidth"
!
      a= d/2.
      b= a/1.732
      c= 2.*b
!
      call plot (x-a,y-b,3)
      call plot (  x,y+c,2)
      call plot (x+a,y-b,2)
      call plot (x-a,y-b,2)
      call plot (x+a,y-b,3)
!
      write(77,*) " 1.0 setlinewidth"
      if(ic.eq.1) return
!
!------------------------
!*  Fill the triangle.
!------------------------
!
      nc=7
      nch= (nc+1)/2
      dx= d/nc
      y2= y -b
!
      do 100 j= 1,nc
      x1= x-a +dx*(j-0.5)
      if(j.le.nch) then
         y1= y +1.732*a*(j-0.5)/nch
      else
         y1= y +1.732*a*(nc-j)/nch 
      end if
!
      call plot (x1,y1,3)
      call plot (x1,y2,2)
  100 continue
!
      call plot (x1,y2,3)
      write(77,*) " 1.0 setlinewidth"
!
      return
      end
!
!
!-------------------------------------------------
       subroutine square (x,y,d,ic)
!-------------------------------------------------
!*  Open square centered at (x,y) /or outer edge.
!
      write(77,*) " 3.0 setlinewidth"
!
      a= d/2.
!
      call plot (x-a,y-a,3)
      call plot (x+a,y-a,2)
      call plot (x+a,y+a,2)
      call plot (x-a,y+a,2)
      call plot (x-a,y-a,2)
      call plot (x-a,y-a,3)
!
      write(77,*) " 1.0 setlinewidth"
      if(ic.eq.1) return
!
!------------------------
!*  Fill the square.
!------------------------
!
      nc=7
      dx= d/nc
      y1= y -a
      y2= y +a
!
      do 100 j= 1,nc
      x1= x-a +dx*(j-0.5)
!
      call plot (x1,y1,3)
      call plot (x1,y2,2)
  100 continue
!
      call plot (x1,y2,3)
      write(77,*) " 1.0 setlinewidth"
!
      return
      end
!
!
!-----------------------------------------------------------------------
      subroutine LPLOT1 (ix,iy,npt1,x,y,ymax,ymin,IL,nxtick,nytick, &
                         LAB1,n1,LAB2,n2,LAB3,n3,ic)
!-----------------------------------------------------------------------
!  <<Warning>>  Order and number of arguments /LPLOT/ have been changed.
!               Also, X (time) is defined for all range.
!               Date: 5/18/96 at MIT.
!***********************************************************************
!   IL=1................ LINEAR PLOT OF (X,Y)
!   IL=2................ LOG10 PLOT OF (X,LOG Y)
!***********************************************************************
      use, intrinsic :: iso_c_binding  ! <-
      implicit none
!
      real(C_float),dimension(npt1) :: x,y,u,v
      real(C_float)  ymax,ymin,XCM(6),YCM(6),PL(6),PR(6),QL(6),QR(6), &
                     xmax,xmin,time,xleng,Edc4,tau_wave
!
      integer(C_int) ix,iy,npt1,IL,nxtick,nytick,n1,n2,n3,ic,iplot
      integer(C_int) npt,isc,i1,j1,i,j,k,ii
!      character(*) isymb  
      character(len=8) lab1,lab2,lab3,label,date*9
!
      real(C_float) HH,dx,dy,scx,scy,x0,x1,x2,x3,x4,y0,y1,y2,y3,y4, &
                    xc,xu,xd,xr,xl,yc,yr,yl,du,uu,smax,smin
!
      COMMON/HEADR1/ label,date
      COMMON/HEADR2/ time,xleng,Edc4,tau_wave
!
      DATA  XCM/12.0, 2*10.00, 3*6.00/,      &
            YCM/15.0,  2*7.50, 3*3.90/,      &
            PL/6.0, 2.0,14.0, 2.0,9.0,16.0/, &
            QL/2.3, 9.5,1.5, 12.9,7.6,2.3/
!
      IPLOT=1
      GO TO 1
!
!-----------------------------------------------------------------------
      entry HPLOT1 (ix,iy,npt1,x,y,ymax,ymin,IL,nxtick,nytick, &
                    LAB1,n1,LAB2,n2,LAB3,n3,ic)
!                                           =1 if shaded histogram.
!-----------------------------------------------------------------------
      IPLOT=2
!
    1 NPT= NPT1
      ISC= 1
!
      do i=1,6
      pr(i)= pl(i) +xcm(i)
      end do
!
      do j= 1,6
      qr(j)= ql(j) +ycm(j)
      end do
!
!                 ******************************************************
!*                **  MAKE A COPY BEFORE THE TOP-LEFT FRAME IS DRAWN. **
!                 ******************************************************
      HH= 0.60
      i1= iabs(ix)
      j1= iabs(iy)
      if(i1.ge.3) go to 10
      if(j1.eq.3 .or. j1.ge.5) go to 10
!                                              ************************
!                                              ** LABEL OF THE PAGE. **
!                                              ************************
!     CALL SYMBOL (0.1,18.7,HH,LABEL,0.,8)
!     CALL SYMBOL (0.1, 0.1,HH,DATE, 0.,9)
      CALL SYMBOL (15.5,0.0,0.7,'t=',0.,2)
      CALL VALUES (16.5,0.0,0.7,time,0.,101)
   10 continue
!
      do i= 1,npt
      u(i)= x(i)
      end do
      xmax= u(npt)
      xmin= u(1)
!                             ************************************
!                             ** THREE-POINT AVERAGE IF IL > 0  **
!                             ************************************
      if(IL.gt.0) then
        v(1)=   y(1)
        v(npt)= y(npt)
        do i= 2,npt-1
        v(i)= 0.33333*(y(i-1)+y(i)+y(i+1))
        end do
      else
        do i= 1,npt
        v(i)= y(i)
        end do
      end if
!                                                *****************
!                                                **  LOG. SCALE **
!                                                *****************
      if(iabs(IL).eq.2) then
         smin= 1.e10
!
         do i= 1,npt
         if(v(i).gt.0.) then
           v(i)= alog10(v(i))
           smin= amin1(smin,v(i))
         else
           v(i)= -1001.
         end if
         end do
!
         do i= 1,npt
         if(v(i).lt.-1000.) v(i)= smin
         end do
!
         call lplmax (v,ymax,ymin,npt)
!        IF(YMAX.GT.0.0) YMAX= YMAX+1.0
      end if
!                                **************************************
!                                ** SET A NEW SCALE AND DRAW A FRAME.**
!                                **************************************
      dx= (xmax-xmin)/xcm(i1)
      dy= (ymax-ymin)/ycm(j1)
      x0= xmin
      y0= ymin
!
      CALL SCALEX (PL(I1),QL(J1),X0,Y0,DX,DY,ISC)
!                                                      *************
!                                                      **  FRAME. **
!                                                      *************
      CALL PLOT (PL(I1),QL(J1),3)
      CALL PLOT (PL(I1),QR(J1),2)
      CALL PLOT (PR(I1),QR(J1),2)
      CALL PLOT (PR(I1),QL(J1),2)
      CALL PLOT (PL(I1),QL(J1),2)
!                                                    ******************
!                                                    **  TICK MARKS. **
!                                                    ******************
      scx= xcm(i1)/(nxtick +1)
      scy= ycm(j1)/(nytick +1)
!
      x0= PL(i1)
      y1= QL(j1)
      y4= QR(j1)
      y2= y1 +0.35
      y3= y4 -0.35
!
      do k= 1,nxtick
      x0= x0 +scx
      CALL PLOT (X0,Y1,3)
      CALL PLOT (X0,Y2,2)
      CALL PLOT (X0,Y3,3)
      CALL PLOT (X0,Y4,2)
      end do
!
      y0= QL(j1)
      x1= PL(i1)
      x4= PR(i1)
      x2= x1 +0.35
      x3= x4 -0.35
!
      do k= 1,nytick
      y0= y0 +scy
      CALL PLOT (X1,Y0,3)
      CALL PLOT (X2,Y0,2)
      CALL PLOT (X3,Y0,3)
      CALL PLOT (X4,Y0,2)
      end do
!                                                     **************
!                                                     ** NUMBERS. **
!                                                     **************
!
!*  gfortran @wat_radial5.f03 
      if(j1.eq.3) then
        CALL NUMBER (PL(I1)-1.3,QL(J1)-0.6,HH,xmin,0.,3)
        CALL NUMBER (PR(I1)-1.3,QL(J1)-0.6,HH,xmax,0.,3)
      end if
!
      CALL NUMBER (PL(I1)-2.70,QL(J1)     ,HH,ymin,0.,101)
      CALL NUMBER (PL(I1)-2.70,QR(J1)-0.30,HH,ymax,0.,101)
!
!                                                     **************
!                                                     **  LABELS. **
!                                                     **************
      xc= 0.5*(PL(i1)+PR(i1))
      xu= xc -0.7
      xd= xc -0.20*N2/2
!
      yr= QR(j1)+0.10
      yl= QL(j1)-1.20
!
      if(ic.eq.0) then
        CALL SYMBOL (XU,YR,0.7,lab1,0.,n1)
        CALL SYMBOL (XD,YL,HH, lab2,0.,n2)
      end if
!
!
      XL= PL(I1)-1.50
      YC= 0.5*(QL(J1)+QR(J1))
      CALL SYMBOL (XL,YC,HH, lab3,0.,n3)
!                                     **********************************
!                                     **  NO PLOT IS MADE IF NPT1 < 0 **
!                                     **********************************
   70 IF(NPT1.LT.0) RETURN
!
      CALL PLOTL (U(1),V(1),ISC,3)
!**
      IF(IPLOT.EQ.1) THEN
         DO I=1,NPT
         CALL PLOTL (U(I),V(I),ISC,2)
         end do
      ELSE
!
         if(IC.ge.0) then
           DO I=1,NPT-1
           CALL PLOTL (U(I+1),V(I)  ,ISC,2)
           CALL PLOTL (U(I+1),V(I+1),ISC,2)
           end do
         else
! Shaded.
           do i= 1,npt-1
           du= u(i+1) -u(i)
!
           do ii= 1,5
           uu= u(i) +du*(ii-1)/5
           CALL PLOTL (UU, 0.0,ISC,3)
           CALL PLOTL (UU,V(I),ISC,2)
           end do 
           end do
         end if
      end if
!**
      CALL PLOTL (U(NPT),V(NPT),ISC,3)
!
      return
      end
!
!
!-----------------------------------------------------------------------
      SUBROUTINE PPLOT (X,Y,NP,ISC)
!-----------------------------------------------------------------------
!
      DIMENSION  X(1),Y(1)
      COMMON/PPLCOM/ NFINE,PL1(10),PR1(10),QL1(10),QR1(10), &
                     XMIN1(10),XMAX1(10),YMIN1(10),YMAX1(10)
!
      IF(NFINE.EQ.0) GO TO 200
      PL= PL1(ISC)
      PR= PR1(ISC)
      QL= QL1(ISC)
      QR= QR1(ISC)
      XMIN= XMIN1(ISC)
      XMAX= XMAX1(ISC)
      YMIN= YMIN1(ISC)
      YMAX= YMAX1(ISC)
!
      RDX= (PR-PL)/(XMAX-XMIN)
      RDY= (QR-QL)/(YMAX-YMIN)
!
      DO 100 I=1,NP
      X0= PL +(X(I)-XMIN)*RDX
      Y0= QL +(Y(I)-YMIN)*RDY
!
      CALL PLOT (X0-0.01,Y0,3)
      CALL PLOT (X0+0.01,Y0,2)
      CALL PLOT (X0,Y0-0.01,3)
      CALL PLOT (X0,Y0+0.01,2)
  100 CONTINUE
!
      RETURN
!
  200 SIZE = 0.01
!
      ISC= 1
      DX= 1./RDX
      DY= 1./RDY
      CALL SCALEX (PL,QL,XMIN,YMIN,DX,DY,ISC)
!
      DO 300 I=1,NP
      CALL PLTPPP(X(I),Y(I),SIZE,ISC)
  300 CONTINUE
!
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      SUBROUTINE PLTPPP(X,Y,R,I)
!-----------------------------------------------------------------------
!
!     ONLY FIRST TIME, CALL PENUP
!
!     (.) PLOT ; R = AMAX1( ABS(XMAX-XMIN),ABS(YMAX-YMIN) ) * 0.001
!
!-----------------------------------------------------------------------
      COMMON/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
!
      CALL PLOT (XL(I)+DXI(I)*(X+R-X0(I)),YL(I)+DYI(I)*(Y+R-Y0(I)),3)
!
      CALL PLOT (XL(I)+DXI(I)*(X+R-X0(I)),YL(I)+DYI(I)*(Y+R-Y0(I)),2)
      CALL PLOT (XL(I)+DXI(I)*(X-R-X0(I)),YL(I)+DYI(I)*(Y+R-Y0(I)),2)
      CALL PLOT (XL(I)+DXI(I)*(X-R-X0(I)),YL(I)+DYI(I)*(Y-R-Y0(I)),2)
      CALL PLOT (XL(I)+DXI(I)*(X+R-X0(I)),YL(I)+DYI(I)*(Y-R-Y0(I)),2)
      CALL PLOT (XL(I)+DXI(I)*(X+R-X0(I)),YL(I)+DYI(I)*(Y+R-Y0(I)),2)
!
      CALL PLOT (XL(I)+DXI(I)*(X+R-X0(I)),YL(I)+DYI(I)*(Y+R-Y0(I)),3)
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      SUBROUTINE PLTQQQ (X,Y,R,I)
!-----------------------------------------------------------------------
!
!     ONLY FIRST TIME, CALL PENUP
!
!     (+) PLOT ; R = AMAX1( ABS(XMAX-XMIN),ABS(YMAX-YMIN) ) * 0.001
!
!-----------------------------------------------------------------------
      COMMON/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
!
      CALL PLOT (XL(I)+DXI(I)*(X+R-X0(I)),YL(I)+DYI(I)*(Y  -Y0(I)),3)
!
      CALL PLOT (XL(I)+DXI(I)*(X+R-X0(I)),YL(I)+DYI(I)*(Y  -Y0(I)),2)
      CALL PLOT (XL(I)+DXI(I)*(X-R-X0(I)),YL(I)+DYI(I)*(Y  -Y0(I)),2)
!
      CALL PLOT (XL(I)+DXI(I)*(X  -X0(I)),YL(I)+DYI(I)*(Y+R-Y0(I)),3)
      CALL PLOT (XL(I)+DXI(I)*(X  -X0(I)),YL(I)+DYI(I)*(Y-R-Y0(I)),2)
!
      CALL PLOT (XL(I)+DXI(I)*(X  -X0(I)),YL(I)+DYI(I)*(Y-R-Y0(I)),3)
!
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      SUBROUTINE SCALEX (XCM,YCM,X00,Y00,DX,DY,ISC)
!-----------------------------------------------------------------------
      COMMON/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
!
      X0(ISC)= X00
      Y0(ISC)= Y00
      DXI(ISC)= 1./DX
      DYI(ISC)= 1./DY
!
      XL(ISC)= XCM
      YL(ISC)= YCM
!
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      SUBROUTINE PLOTL (X,Y,ISC,IPL)
!-----------------------------------------------------------------------
      COMMON/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
!
      XCM= XL(ISC) +DXI(ISC)*(X -X0(ISC))
      YCM= YL(ISC) +DYI(ISC)*(Y -Y0(ISC))
!
      CALL PLOT (XCM,YCM,IPL)
!
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      SUBROUTINE PLOTB (X,Y,ISC,IPL)
!-----------------------------------------------------------------------
      COMMON/GSCALE/ X0(10),Y0(10),XL(10),YL(10),DXI(10),DYI(10)
!
      XCM= XL(ISC) +DXI(ISC)*(X -X0(ISC))
      YCM= YL(ISC) +DYI(ISC)*(Y -Y0(ISC))
!
      call plot (xcm,ycm,ipl)
!
!     CALL DASHP (XCM,YCM,IPL)
!
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      SUBROUTINE VALUES(X,Y,HEIGHT,VAL,THETA,IFMAT)
!-----------------------------------------------------------------------
!  << VALUES >>
!     1. FUNCTION
!        (1) TO DRAW VARIABLE
!     2. ARGUMENTS   (SIZE)   (I/O)     (MEANING)
!        (1) X,Y               (I)       ABSOLUTE COORDINATE VALUE
!        (2) HEIGHT            (I)       DRAW OUT SIZE ON PAPER
!        (3) VAL               (I)       VARIABLE
!        (4) THETA             (I)       ANGLE
!        (5) IFMAT             (I)       FORMAT TYPE
!     3. CALLED BY
!             (** NOTHING **)
!     4. CALLS
!             (** NUMBER **)
!             (** SYMBOL **)
!-----------------------------------------------------------------------
!        IFMAT = (N100)*100 + KETA
!        N100 = 0 : INTEGER FORMAT
!        N100 = 1 : F FORMAT ::  NUMBER(X,Y,HEIGHT,VAL,THETA,KETA)
!        N100 = 2 : E FORMAT ::
!        N100 = 3 : POWER OF TEN FORMAT
!        N100 = OTHEWISE : NOT WRITE OUT
!-----------------------------------------------------------------------
!
      REAL*4 VAL
      CHARACTER CHR13*13,CHR12*12,CHR3*3
      CHARACTER*1 MINUS,ZERO,BLANK
      PARAMETER(RATIO = 6./7. )
      DATA MINUS/'-'/,ZERO/'0'/,BLANK/' '/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF (IFMAT.LT.0) RETURN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      N100 = IFMAT/100
      KETA = IFMAT - N100*100
!
      IF (N100.EQ.0) THEN
        CALL NUMBER(X,Y,HEIGHT,VAL,THETA,-1)
      ELSE IF (N100.EQ.1) THEN
        CALL NUMBER(X,Y,HEIGHT,VAL,THETA,KETA)
      ELSE IF (N100.EQ.2) THEN
        CHR13 = '             '
        CHR12 = '            '
        IF (KETA.EQ.0) THEN
          WRITE(CHR13,'(1PE13.6)') VAL
          CHR12(1:4) = CHR13(1:3)//'E'
          NUMSYM = 4
        ELSE
          KETA = KETA + 1
          IF (VAL.LT.0.) THEN
            CHRVAL = VAL - 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+3) = CHR13(1:KETA+2)//'E'
            NUMSYM = KETA + 3
          ELSE IF (VAL.EQ.0) THEN
            CHRVAL = VAL
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+3) = CHR13(1:KETA+2)//'E'
            NUMSYM = KETA + 3
          ELSE
            CHRVAL = VAL + 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+2) = CHR13(2:KETA+2)//'E'
            NUMSYM = KETA + 2
          END IF
        END IF
        CHR3 = '   '
!
        IF (CHR13(11:11) .EQ. MINUS) THEN
          IF (CHR13(12:12) .EQ. ZERO  .OR. &
              CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:2) = '-'//CHR13(13:13)
          ELSE
            CHR3(1:3) = '-'//CHR13(12:13)
          END IF
          NUMSY1 = 3
        ELSE
          IF (CHR13(12:12) .EQ. ZERO  .OR. &
              CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:1) = CHR13(13:13)
            NUMSY1 = 1
          ELSE
            CHR3(1:2) = CHR13(12:13)
            NUMSY1 = 2
          END IF
        END IF
        AKAKU = 2. * 3.1415927 / 360.
        COST = COS(THETA*AKAKU)
        CALL SYMBOL(X,Y,HEIGHT,CHR12,THETA,NUMSYM)
        CALL SYMBOL(999.,999.,HEIGHT,CHR3,THETA,NUMSY1)
      ELSE IF (N100.EQ.3) THEN
        CHR13 = '             '
        CHR12 = '            '
        IF (KETA.EQ.0) THEN
          WRITE(CHR13,'(1PE13.6)') VAL
          CHR12(1:6) = CHR13(1:3)//'X10'
          NUMSYM = 6
        ELSE
          KETA = KETA + 1
          IF (VAL.LT.0.) THEN
            CHRVAL = VAL - 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+5) = CHR13(1:KETA+2)//'X10'
            NUMSYM = KETA + 5
          ELSE
            CHRVAL = VAL + 5.*10**FLOAT(-KETA)
            WRITE(CHR13,'(1PE13.6)') CHRVAL
            CHR12(1:KETA+4) = CHR13(2:KETA+2)//'X10'
            NUMSYM = KETA + 4
          END IF
        END IF
        CHR3 = '   '
!
        IF (CHR13(11:11) .EQ. MINUS) THEN
          IF (CHR13(12:12) .EQ. ZERO  .OR. &
              CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:2) = '-'//CHR13(13:13)
          ELSE
            CHR3(1:3) = '-'//CHR13(12:13)
          END IF
          NUMSY1 = 3
        ELSE
          IF (CHR13(12:12) .EQ. ZERO  .OR. &
              CHR13(12:12) .EQ. BLANK) THEN
            CHR3(1:1) = CHR13(13:13)
            NUMSY1 = 1
          ELSE
            CHR3(1:2) = CHR13(12:13)
            NUMSY1 = 2
          END IF
        END IF
        AKAKU = 2. * 3.1415927 / 360.
        COST = COS(THETA*AKAKU)
        SINT = SIN(THETA*AKAKU)
        CALL SYMBOL(X,Y,HEIGHT,CHR12,THETA,NUMSYM)
!
!                                             *******************
!                                             ** EXPONENT PART **
!                                             *******************
!
        H2 = HEIGHT * 5./7.
        X1 = (NUMSYM+1)* HEIGHT * RATIO
        Y1 = HEIGHT * 4./7.
        IF (ABS(THETA).LT.1E-04) THEN
          X1 = X + X1
          Y1 = Y + Y1
        ELSE
          X2 =     X1 * COST - Y1 * SINT
          Y1 = Y + X1 * SINT + Y1 * COST + H2*COST
          X1 = X + X2                    - H2*SINT
        END IF
        CALL SYMBOL(X1,Y1,H2,CHR3,THETA,NUMSY1)
      END IF
      RETURN
      END
!
!
!-----------------------------------------------------------------------
      SUBROUTINE CIRC(XC,YC,RLENG,NT)
!-----------------------------------------------------------------------
!
      PI2 = 2.*3.1415927
      DT0 = PI2/FLOAT(NT)
      CALL PLOT(XC+RLENG,YC,3)
      DO 1000 ICIRC = 1,NT
        XX = RLENG*COS(DT0*ICIRC) + XC
        YY = RLENG*SIN(DT0*ICIRC) + YC
        CALL PLOT (XX,YY,2)
 1000 CONTINUE
      RETURN
      END
!
!
!***************************************************************
!*     This program package generates a UNIX postscript        *
!*     graphic file when called by calcomp-compatible          *
!*     /plot23.f/.                                             *
!***************************************************************
!----------------------------------------------------------
!      PostScript header by fortran
!        T. Ogino (Nagoya University) February 27, 1992
!      Modified to conform GSIPP commands
!        Motohiko Tanaka (NIFS)       November 23, 1993
!
!----------------------------------------------- 5/27/96 -------
!     This PS-Adobe-2.0 header allows us full paging features in
!     the Ghostview.  To scroll up the page (backward), click the 
!     page number and press two buttons of mouse simultaneously.
!
!     Consult: A.Saitou (Kyoto U.)  The definition of /@eop  
!    needs stroke for line drawings (not in the TeX header).
!---------------------------------------------------------------
       subroutine gopen (nframe)
!----------------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
       common/pages/  ipage,nfrm
!
!*  This is an Adobe-2.0 postscript file.
!
       write(77,10)
   10  format('%!PS-Adobe-2.0',/      &
              '%%Pages: (atend)',/    &
              '%%PageOrder: Ascend',/ &
              '%%EndComments',/       &
              '%%BeginDocument')
!
!%%%%%%%%%%%%%%%%%%% Procedure Defintions %%%%%%%%%%%%%%%%%%%%%%%%%%
!
      write(77,11) 
   11 format('%%BoundingBox: 150. 400. 550. 600.')
!
      write(77,21) 
   21 format('/l {lineto} bind def  % x y l -- line to position',/ &
             '/m {moveto} bind def  % x y m -- move to position')
!
      write(77,23) 
   23 format('/tr {/Times-Roman findfont} bind def',/ &
             '/sf {scalefont} bind def',/ &
             '/se {setfont} bind def',/   &
             '/ro {rotate}  bind def',/   &
             '/tl {translate} bind def',/ &
             '/sc {scale} bind def')
!
      write(77,24) 
   24 format('/@bop          % @bop -- begin the a new page',/ &
             '{erasepage newpath initgraphics',/ &
             '/SaveImage save def',/             &
             '} bind def')
!
      write(77,25) 
   25 format('/@eop          % @eop -- end a page',/ &
             '{stroke showpage',/   &
             ' SaveImage restore',/ &
             '} bind def')
!
      write(77,26) 
   26 format('/@end          % @end -- done the whole shebang',/ &
             ' /end load def')
!
      write(77,27) 
   27 format('/dir 0 def')
!
      write(77,29) 
   29 format('/s             % string s -- show the string',/ &
             '{dir 1 eq',/ &
             ' {gsave currentpoint translate 90 rotate 0 0 moveto',/ &
             ' show grestore}',/ &
             ' {show} ifelse',/  &
             '} bind def')
!
      write(77,31)
   31 format('%%EndDocument',/  &
             '%%EndProlog',/    &
             '%%BeginSetup',/   &
             '/Resolution 300 def',/ &
             '/#copies 1 def',/ &
             '%%EndSetup')
!
!%%%%%%%%%%%%%%%%%%% End of the header %%%%%%%%%%%%%%%%%%%%%%%%%%
!
!*  initiate the page one.
!
       nfrm = nframe
!
       ipage = 1
       write(77,12) ipage,ipage
   12  format('%%Page:',1x,i2,1x,i2)
!
       write(77,30) 
   30  format('%%BeginPageSetup',/ &
              '%%EndPageSetup',/   &
              '@bop')
!
!
!*  Set magnifying factor (GSIPP to Sun coordinate).
!   Rotate and translate to output on A4-L paper.
!      Left corner ...... (  0.,  0.)
!      Right corner ..... (600.,780.)
!
       xcm=  25.
       xwc= 700.
       fmag= xwc/xcm
!
       write(77,*) '90.0 ro'
       write(77,*) '50.0 -550.0 tl'
!
!*  If nfrm=4, four frames in a page (top-left frame).
!
       if(nfrm.eq.1) then
          write(77,*) '1.00 1.00 sc'
       else
          write(77,*) '0.50 0.50 sc'
          write(77,*) '0.0 550.0 tl'
       end if
!
       return
       end
!
!
!-----------------------------
       subroutine gclose
!-----------------------------
       call plote
       return
       end
!
!
!-----------------------------
       subroutine plote
!-----------------------------
       write(77,10) 
   10  format('@eop')
       return
       end
!
!
!-----------------------------------------
       subroutine chart
!-----------------------------------------
!*     Four frames in a page (if nfrm=4).
       common/pages/ ipage,nfrm
!
!
       ipage = ipage +1
       loc= mod(ipage-1,nfrm)
!
!*  Frame 1: open a new page.
!
       if(loc.eq.0) then
          call plote
!
          if(nfrm.eq.1) lpage= ipage
          if(nfrm.ne.1) lpage= (ipage+3)/4
!
          write(77,10) 
   10     format('%%PageTrailer    % Need for the page count')
!
          write(77,20) lpage,lpage
   20     format('%%Page:',1x,i2,1x,i2)
!
          write(77,30) 
   30     format('%%BeginPageSetup',/ &
                 '%%EndPageSetup',/   &
                 '@bop')
!
          write(77,*) '90.0 ro'
          write(77,*) '50.0 -550.0 tl'
!
          if(nfrm.eq.1) then
             write(77,*) '1.00 1.00 sc'
          else
             write(77,*) '0.50 0.50 sc'
             write(77,*) '0.0  550.0 tl'
          end if
!
          return
       end if
!
!
!-----------------------------------------------------
!      First cancel the previous translation, then
!      make a new translation (scale factor alive).
!-----------------------------------------------------
!*   Frames 2-4:
!
       if(loc.eq.1) then
          write(77,*) '  0.0 -550.0 tl'
          write(77,*) '700.0  550.0 tl'
       end if
!
       if(loc.eq.2) then
          write(77,*) '-700.0 -550.0 tl'
          write(77,*) '   0.0    0.0 tl'
       end if
!
       if(loc.eq.3) then
          write(77,*) '  0.0 0.0 tl'
          write(77,*) '700.0 0.0 tl'
       end if
!
       return
       end
!
!
!------------------------------------
       subroutine factor(fct)
!------------------------------------
       write(77,10) fct,fct
   10  format(f6.2,1x,f6.2,' sc')
       return
       end
!
!
!------------------------------------
       subroutine newpen (ip)
!------------------------------------
       i1=(ip-1)/2
       i2=ip-2*i1
       write(77,*) 'sn'
       pi1=0.40*float(i1-1)
       write(77,30) pi1
   30  format(f3.1,' sl')
       if(i2.ne.1) then
       write(77,*) '[2 2] 0 sd'
       endif
       return
       end
!
!
!-----------------------------
       subroutine linee
!-----------------------------
       write(77,*) 'st'
       return
       end
!
!
!------------------------------------
       subroutine plot (x0,y0,ip)
!------------------------------------
       use, intrinsic :: iso_c_binding  ! <-
       implicit none
!
       real(C_float) x0,y0,x,y,h
       integer(C_int) ip,n
!
       x= x0
       y= y0
       h= 0.
       n= 777
       call sunscl (x,y,h,n)
!
       if(ip.eq.3)  write(77,10) x,y
       if(ip.eq.2)  write(77,20) x,y
       if(ip.eq.-3) write(77,30) x,y
       if(ip.eq.-2) write(77,40) x,y,x,y
   10  format(f5.1,1x,f5.1,' m')
   20  format(f5.1,1x,f5.1,' l')
   30  format(f5.1,1x,f5.1,' tl')
   40  format(f5.1,1x,f5.1,' l sn',1x,f5.1,1x,f5.1,' tl')
!       write(77,*) 'st'
!
       return
       end
!
!
!-------------------------------------------------
       subroutine symbol (x0,y0,h0,isymb,ang,n0)
!-------------------------------------------------
       use, intrinsic :: iso_c_binding  ! <-
       implicit none
!
       integer(C_int) iffig
       common/comfig/ iffig
!
       real(C_float) x0,y0,x,y,h
       real(C_float) h0,ang
       character(*) isymb   !!! nec: hitsuyo
       integer(C_INT) n0,n,i
!
       character    ica*80,ich(80)*1
       equivalence (ica,ich(1))
!
       x= x0
       y= y0
       h= h0
       n= n0
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!*
       ica= isymb
       write(77,*) '(',(ich(i),i=1,n),') s'
!
       return
       end
!
!
!-----------------------------------------------
       subroutine number (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       use, intrinsic :: iso_c_binding  ! <-
       implicit none
!
       integer(C_int) iffig
       common/comfig/ iffig
!
       real(C_float) x0,y0,h0,anu,ang,x,y,h
       integer(C_INT) n0,n
       character  isymb*9
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
!
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!
       if(abs(anu).gt.1.e+1 .or.  &
          abs(anu).lt.1.e-1) then
         write(isymb,31) anu
   31    format(1pe9.2)
       else
         write(isymb,32) anu
   32    format(f7.2)
       end if
!
       if(.true.) go to 300
       if(abs(anu).lt.10000.) then  ! 5 digits
         if(abs(anu).gt.0.1) then
           write(isymb,40) anu
   40      format(f6.1)
         else
           if(abs(anu).gt.0.001) then  ! f6.3
             write(isymb,41) anu
   41        format(f6.3)
           else
             if(abs(anu).gt.0.001) then  ! f6.3
               write(isymb,42) anu   ! 1pe9.2
   42          format(1pe9.2)
             else
               write(isymb,40) anu   ! f6.1
             end if
           end if
         end if
!
       else
         if(abs(anu).lt.100000.) then
           write(isymb,51) anu     ! f7.1
   51      format(f7.1)
         else
           write(isymb,52) anu     ! e9.2
   52      format(1pe9.2)
         end if
       end if
  300  continue
!
       write(77,*) '(',isymb,') s'
!
       return
       end
!
!
!-----------------------------------------------
       subroutine number2 (x0,y0,h0,anu,ang,n0)
!-----------------------------------------------
       use, intrinsic :: iso_c_binding  ! <-
       implicit none
!
       integer(C_int) iffig
       common/comfig/ iffig
!
       real(C_float) x0,y0,h0,anu,ang,x,y,h
       integer(C_INT) n0,n
       character  isymb*6
!
       x= x0
       y= y0
       h= h0
       n= 777
       call sunscl (x,y,h,n)
!
       write(77,*) 'tr'
       write(77,10) h
   10  format(f5.1,' sf')
       write(77,*) 'se'
!
       write(77,20) x,y
   20  format(f5.1,1x,f5.1,' m')
       write(77,30) ang
   30  format(f5.1,' ro')
!
       if(n0.eq.1) write(isymb,41) anu
       if(n0.eq.2) write(isymb,42) anu
   41  format(f6.1)
   42  format(f6.2)
!
       write(77,*) '(',isymb,') s'
!
       return
       end
!
!
!---------------------------------------------------
       subroutine sunscl (x,y,h,n)
!---------------------------------------------------
       common/convsn/ fmag,x0,y0,h0,n0
!
       if(x.eq.999.) then
         x= x0 +iabs(n0)*h0
       else
         x= fmag*x
         x0= x
       end if
!
       if(y.eq.999.) then
         y= y0
       else
         y= fmag*y
         y0= y
       end if
!
       h= fmag*h
       h0= h
       if(n.ne.777) n0= n
!
       return
       end
