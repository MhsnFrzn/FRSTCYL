c      SUBROUTINE initial(dt0,ch_fin,m_fin,y_toto,z_toto)
      SUBROUTINE initial(dt0,ch_fin,m_fin)
#include "common.inc"
      REAL*8 dt0
cccccccccccccccccccccccccccccc GA cccccccccccccccc
      REAL*8 aa,bb,cc,dd,ee,ff,a,b,c,d,e,f,alfa,beta,nznz
      REAL*8 xxk,yyk,c1000,c2000
      INTEGER lwork,i22,i33,ipiv1(nz)
      PARAMETER (lwork=nz*(nz+4))
      REAL*8 wavex_ratio(nx),wavey_ratio(ny),work(lwork)
      REAL*8 R1z(nz,nz),L1z(nz,nz),aalpha
      INTEGER css,scs,ssc
      REAL*8 D2z(nz,nz),cn_drdr_dum(nz,nz),e_drdr_dum(nz,nz)
      REAL*8 ci_drdr_dum(nz,nz),cn_drdr_inv(nz,nz),D2z2(nz,nz)
      CHARACTER*32 anum
cccccccccccccccccccccccccccccc GA 
      INTEGER ch_fin,m_fin
c
c     set up the matrices and ffts, initialize the field variables
c
      INTEGER i,j,k,pt,n,l
      REAL*8 pin,hpin,sum,c100,c200,xk,yk,d00,d0n,dn0,dnn,iden
     &      ,inzm,tinzm,hinzm
      REAL*8 condo,condi,emax,emin
      REAL*8 aux_rr,h0,sumtemp !,theta
      REAL*8 thickness,center_tt,hwidth_tt
      
cccccccccccccccccccccccccccc GA cccccccccccccccc
c 1 to turn on ibm, 0 for off
c      flagibm=1 
c 1 for compact 0 for spectral in x direction
      css=0
c 1 for compact 0 for spectral in y direction
      scs=0
c 1 for compact 0 for spectral in z direction
      ssc=0
      open(unit=100, file='schemes')
cccccccccccccccccccccccccccc GA
      if (FLTIMER) call begin_timer(10)
c  set write and read format
#ifdef IFC
C      open(21,file=name,form='binary')
      format_mode='binary'
#else
C      open(21,file=name,form='unformatted')
      format_mode='unformatted'
#endif
c
c  set some variables
      pi = 4.d0*atan(1.d0)
      ipin = dble(nzm)/pi
c
c  rlenxm rleny and rlenz are input on chnl.anz_ini file
c  here recompute using pi for horizontal box size for 
c  channel flow 
      rlenx = rlenx   !4.d0*pi/3.d0
      rleny = rleny  !4.d0*pi
      write(26,*)'Internally computed rlenx and rleny'
c
      deltax = rlenx/dble(nx)
      deltay = rleny/dble(ny)
      idx = 1.d0/deltax
      idy = 1.d0/deltay
      dt0 = 0.d0
cccccccccccccccccccccccccccccccccccccc GA cccccccccccccccc
      do i=1,nx0
         xets(i)=(i-1)*deltax
      enddo
      do j=1,ny0
         yets(j)=(j-1)*deltay
      enddo
cccccccccccccccccccccccccccccccccccccc GA
c
c  checking some values for mean computations
      if(np.eq.1.and.dnp.ne.1)then
         write(26,*)'Wrong value for dnp:',dnp
         write(26,*)'if np=1 => dnp=1'
         write( *,*)'Wrong value for dnp:',dnp
         write( *,*)'if np=1 => dnp=1'
         stop
      endif
c
c  compute dt internally for oscillatory flow
      if(dt.lt.0.d0)then 
         dt=2.d0*pi/dble(np)/dble(dnp)
         write(26,*)'computing dt internally: dt = ',dt
      endif
c 
c  make sure ws is magnitude of settling velocity
      ws=abs(ws)
c
cNOTE: see below for gravity computation
c
c     zero multidimensional arrays
c---------------------------------------------------------------------
c
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP DO PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny0
          do i=1,nx0
            u(  i,j,k) = 0.d0
            v(  i,j,k) = 0.d0
            w(  i,j,k) = 0.d0
            p(  i,j,k) = 0.d0
cccccccccccccccccccccccccccccccccccccc GA cccccccccccccccc
            dpx(i,j,k) = 0.d0
            dpy(i,j,k) = 0.d0
            dpz(i,j,k) = 0.d0
cccccccccccccccccccccccccccccccccccccc GA
            ru( i,j,k) = 0.d0
            rv( i,j,k) = 0.d0
            rw( i,j,k) = 0.d0
            rp( i,j,k) = 0.d0
            h1( i,j,k) = 0.d0
            h2( i,j,k) = 0.d0
            h3( i,j,k) = 0.d0
            rtt(i,j,k) = 0.d0
            ru2( i,j,k) = 0.d0
            rv2( i,j,k) = 0.d0
            rw2( i,j,k) = 0.d0
            rp2( i,j,k) = 0.d0
            rtt2(i,j,k) = 0.d0
          enddo
        enddo
      enddo
C$OMP END DO NOWAIT
      if (FLTHRM) then
C$OMP DO PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny0
               do i=1,nx0
                  tt( i,j,k) = 0.d0
                  htt(i,j,k) = 0.d0
               enddo
            enddo
         enddo
C$OMP END DO NOWAIT
      endif
c
      if (FLTHRM2) then
C$OMP DO PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny0
               do i=1,nx0
                  tt2( i,j,k) = 0.d0
                  htt2(i,j,k) = 0.d0
               enddo
            enddo
         enddo
C$OMP END DO NOWAIT
      endif
C$OMP SECTIONS PRIVATE(i,j)
      do j=1,nyh
         do i=1,nx0
            pbxnb(i,j) = 0.d0
            pcxnb(i,j) = 0.d0
            pbynb(i,j) = 0.d0
            pcynb(i,j) = 0.d0
            pbxob(i,j) = 0.d0
            pcxob(i,j) = 0.d0
            pbyob(i,j) = 0.d0
            pcyob(i,j) = 0.d0
            pbxnt(i,j) = 0.d0
            pcxnt(i,j) = 0.d0
            pbynt(i,j) = 0.d0
            pcynt(i,j) = 0.d0
            pbxot(i,j) = 0.d0
            pcxot(i,j) = 0.d0
            pbyot(i,j) = 0.d0
            pcyot(i,j) = 0.d0
         enddo
      enddo
c
      do j=1,ny0
         do i=1,nx0
            ubct(   i,j)=0.d0   !=ugt if non-homog, see uvci_mimi
            ubcb(   i,j)=0.d0   !=ugt if non-homog, see uvci_mimi
            vbct(   i,j)=0.d0   !=vgt if non-homog, see uvci_mimi
            vbcb(   i,j)=0.d0   !=vgt if non-homog, see uvci_mimi
            wbct(   i,j)=0.d0   !=wgt if non-homog, see  wci_mimi
            wbcb(   i,j)=0.d0   !=wgt if non-homog, see  wci_mimi
            ttbct(  i,j)=ttgt(1)
            ttbcb(  i,j)=ttgb(1)
            tt2bct(  i,j)=tt2gt(1)
            tt2bcb(  i,j)=tt2gb(1)
            deposit(i,j)=0.d0
            erosion(i,j)=0.d0
         enddo
      enddo
c
      do j=1,nz
         do i=1,nz0
            chb(    i,j) = 0.d0
            chbd1(  i,j) = 0.d0
            chbd2(  i,j) = 0.d0
            cn_drdr(i,j) = 0.d0
            ci_drdr(i,j) = 0.d0
            cn_nodr(i,j) = 0.d0
            ci_nodr(i,j) = 0.d0
         enddo
         do i=1,nyl0
            gy(i,j) = 0.d0
            gz(i,j) = 0.d0
         enddo
      enddo
      do j=1,nzmm
         do i=1,nzmm0
            cn_nono(i,j)   = 0.d0
            ci_nono(i,j)   = 0.d0
            uvcn_mimi(i,j) = 0.d0
            uvci_mimi(i,j) = 0.d0
            wcn_mimi(i,j)  = 0.d0
            wci_mimi(i,j)  = 0.d0
            ttcn_mimi(i,j) = 0.d0
            ttci_mimi(i,j) = 0.d0
            tt2cn_mimi(i,j) = 0.d0
            tt2ci_mimi(i,j) = 0.d0
         enddo
      enddo
C$OMP SECTION
      do i=1,np
         do k=1,nz
            um(   k,i) = 0.d0
            vm(   k,i) = 0.d0
            wm(   k,i) = 0.d0
            ttm(  k,i) = 0.d0
            u2m(  k,i) = 0.d0
            v2m(  k,i) = 0.d0
            w2m(  k,i) = 0.d0
            tt2m( k,i) = 0.d0
            u3m(  k,i) = 0.d0
            v3m(  k,i) = 0.d0
            w3m(  k,i) = 0.d0
            tt3m( k,i) = 0.d0
            u4m(  k,i) = 0.d0
            v4m(  k,i) = 0.d0
            w4m(  k,i) = 0.d0
            tt4m( k,i) = 0.d0
            uvm(  k,i) = 0.d0
            uwm(  k,i) = 0.d0
            uttm( k,i) = 0.d0
            vwm(  k,i) = 0.d0
            vttm( k,i) = 0.d0
            wttm( k,i) = 0.d0
c     
            uxm(  k,i) = 0.d0
            uym(  k,i) = 0.d0
            uzm(  k,i) = 0.d0
            vxm(  k,i) = 0.d0
            vym(  k,i) = 0.d0
            vzm(  k,i) = 0.d0
            ttxm( k,i) = 0.d0
            ttym( k,i) = 0.d0
            ttzm( k,i) = 0.d0
c     
            ux2m( k,i) = 0.d0
            uy2m( k,i) = 0.d0
            uz2m( k,i) = 0.d0
            vx2m( k,i) = 0.d0
            vy2m( k,i) = 0.d0
            vz2m( k,i) = 0.d0
            ttx2m(k,i) = 0.d0
            tty2m(k,i) = 0.d0
            ttz2m(k,i) = 0.d0
c     
            ux3m( k,i) = 0.d0
            uy3m( k,i) = 0.d0
            uz3m( k,i) = 0.d0
            vx3m( k,i) = 0.d0
            vy3m( k,i) = 0.d0
            vz3m( k,i) = 0.d0
            ttx3m(k,i) = 0.d0
            tty3m(k,i) = 0.d0
            ttz3m(k,i) = 0.d0
c     
            ux4m( k,i) = 0.d0
            uy4m( k,i) = 0.d0
            uz4m( k,i) = 0.d0
            vx4m( k,i) = 0.d0
            vy4m( k,i) = 0.d0
            vz4m( k,i) = 0.d0
            ttx4m(k,i) = 0.d0
            tty4m(k,i) = 0.d0
            ttz4m(k,i) = 0.d0
            pwm(  k,i) = 0.d0
            pm(   k,i) = 0.d0
            u2wm( k,i) = 0.d0
            v2wm( k,i) = 0.d0
            uuzm( k,i) = 0.d0
            uwxm( k,i) = 0.d0
            vvzm( k,i) = 0.d0
            vwym( k,i) = 0.d0
            wwzm( k,i) = 0.d0
         enddo
      enddo
C$OMP SECTION
      do i=1,np
         do k=1,nz
            do j=1,nxh
               uespxm(j,k,i)=0.d0
               vespxm(j,k,i)=0.d0
               wespxm(j,k,i)=0.d0
c               udspxm(j,k,i)=0.d0
c               vdspxm(j,k,i)=0.d0
c               wdspxm(j,k,i)=0.d0
            enddo
            do j=1,nyh
               uespym(j,k,i)=0.d0
               vespym(j,k,i)=0.d0
               wespym(j,k,i)=0.d0
c               udspym(j,k,i)=0.d0
c               vdspym(j,k,i)=0.d0
c               wdspym(j,k,i)=0.d0
            enddo
         enddo
      enddo
C$OMP END SECTIONS
C$OMP END PARALLEL
c
c     define the coefficients for the r-k time stepping
c
c     coefficient of old advection term
c
      cfnl1v(1) = 0.d0
      cfnl1v(2) = -5.d0/9.d0
      cfnl1v(3) = -153.d0/128.d0
c
c     coefficient of advection term
c
      cfnl2v(1) = 1.d0/3.d0
      cfnl2v(2) = 15.d0/16.d0
      cfnl2v(3) = 8.d0/15.d0
c
c     coefficient of diffusion term
c
      cfdifv(1) = dt/6.d0
      cfdifv(2) = dt*5.d0/24.d0
      cfdifv(3) = dt/8.d0
c
c     sub-step sizes for RK3 stepping
c
      cfdelv(1) = 1.d0/3.d0
      cfdelv(2) = 5.d0/12.d0
      cfdelv(3) = 1.d0/4.d0
c
c     accumulation of sub-step sizes
c
      cfcumv(1) = 1.d0/3.d0
      cfcumv(2) = 3.d0/4.d0
      cfcumv(3) = 1.d0
c
c     set erosion coefficients
      cferos1=1.d0
      cferos2=0.6d0
      if(rey_p.le.2.34)then
         cferos1=0.586d0
         cferos2=1.23d0
      endif
c
c     set up fourier apparatus
c     initialize fft's
      call initfft
c
      c100 = 2.d0*pi/rlenx
      c200 = 2.d0*pi/rleny
c
c     wavex is the wavenumber, for a complex DFT
c     wavexs is the square of the wavenumber
c     the subarray 1:nxhp is used for a real DFT
c
      do i=1,nx
        wavex(i) = 0.d0
        wavexs(i) = 0.d0
      enddo
cccccccccccccccccccccccccccccccccccccccccccc GA cccccccccccccccccccccc
c     wavex_ratio is the ratio that can change spectral to compact in x direction
c     css=1 will have compact scheme in the x direction
      if(css.eq.1)then
        do i=1,nxh
          wavex_ratio(i)=0.d0
        enddo
c
        c1000=2.d0*pi/dble(nx)
        do i=2,nxh
          xk = c100*dble(i-1)
          xxk = c1000*dble(i-1)
          wavex_ratio(i)=1.d0/xxk*(14.d0/9.d0*sin(xxk)
     &     +1.d0/18.d0*sin(2.d0*xxk))/(1.d0+2.d0/3.d0*cos(xxk))
          wavex(i) = xk*wavex_ratio(i)
          wavexs(i) = wavex(i)**2
          wavex(nx+2-i) = -xk*wavex_ratio(i)
          wavexs(nx+2-i) = wavex(nx+2-i)**2
        enddo
      write(100,*) 'x compact'
      endif
ccccccccccccccccccccccccccccccccccccccccccccc GA       
c     css=0 will have spectral scheme in the x direction
      if(css.eq.0)then
        do i=2,nxh
          xk = c100*dble(i-1)
          wavex(i)=xk
          wavexs(i) = xk**2
          wavex(nx+2-i) = -xk
          wavexs(nx+2-i) = xk**2
        enddo
      write(100,*) 'x spectral'
      endif
c
      do j=1,ny
        wavey(j) = 0.d0
        waveys(j) = 0.d0
      enddo
ccccccccccccccccccccccccccccccccccccccccccccc GA ccccccccccccccccccccc
c     wavey_ratio is the ratio that can change spectral to compact in y direction
c     scs=1 will have compact scheme in the y direction
      if(scs.eq.1)then
        do j=1,nyh
          wavey_ratio(j)=0.d0
        enddo
c
        c2000=2.d0*pi/dble(ny)
        do j=2,nyh
          yk = c200*dble(j-1)
          yyk = c2000*dble(j-1)
          wavey_ratio(j)=1.d0/yyk*(14.d0/9.d0*sin(yyk)
     &      +1.d0/18.d0*sin(2.d0*yyk))/(1.d0+2.d0/3.d0*cos(yyk))
          wavey(j) = yk*wavey_ratio(j)
          waveys(j) = wavey(j)**2
          wavey(ny+2-j) = -yk*wavey_ratio(j)
          waveys(ny+2-j) = wavey(ny+2-j)**2
        enddo
      write(100,*) 'y compact'
      endif
ccccccccccccccccccccccccccccccccccccccccccccc GA
c     scs=0 will have spectral scheme in the y direction
      if(scs.eq.0)then
        do j=2,nyh
          yk = c200*dble(j-1)
          wavey(j) = yk
          waveys(j) = yk**2
          wavey(ny+2-j) = -yk
          waveys(ny+2-j) = yk**2
        enddo
      write(100,*) 'y spectral'
      endif
c
      write(26,107)(i,wavex(i),i=1,nx)
      write(26,108)(i,wavey(i),i=1,ny)
 107  FORMAT('   wavex  ',/,129(2x,i4,2x,1pe15.7,/))
 108  FORMAT('   wavey  ',/,129(2x,i4,2x,1pe15.7,/))
c
c     compute vertical grid points
cNOTE: here the Gauss-Lobatto points are computed for the vertical
c      direction and are in the interval [-1,1]
c      See below for stretching in the interval [-rlenz/2,rlenz/2]
c
      pin  = pi/dble(nzm)
      hpin = 0.5d0*pin
ccccccccccccccccccccccccccccccccccccccccccccc GA ccccccccccccccccccccc
c set aalpha to 0 for GL pts. Alpha aproaches 1 for uniform. Cannot be 1
      aalpha=0.d0
c      aalpha=0.99d0
      do k = 1,nz/2
        if(aalpha.eq.0)then
           zets(k) = cos(dble(k-1)*pin)
           zets(nz+1-k) = -zets(k)
        else
           zets(k)=asin(aalpha*cos(dble(k-1)*pin))/asin(aalpha)
           zets(nz+1-k) = -zets(k)
        endif
ccccccccccccccccccccccccccccccccccccccccccccc GA
        sum = 1.d0
        if (mod(nz,2).eq.1) sum = 1.d0 + (-1)**k / dble(nzm**2-1)
        do i = 1,nz/2-1
          sum = sum - 2.d0*cos(2.d0*pin*dble(i*(k-1)))/dble(4*i**2-1)
        enddo
        sum = sum/dble(nzm)
        if (k.eq.1) sum = 0.5d0*sum
        wzet(k) = sum
        wzet(nz+1-k) = sum
      enddo
      if (mod(nz,2).eq.1) then
        k = nz/2 + 1
        zets(k) = 0.d0
        sum = 1.d0 + (-1)**k / dble(nzm**2-1)
        do i = 1,nz/2-1
          sum = sum - 2.d0*(-1)**i / dble(4*i**2-1)
        enddo
        wzet(k) = sum/dble(nzm)
      endif
      do k = 1,nzm/2
ccccccccccccccccccccccccccccccccccccccccccccc GA ccccccccccccccccccccc
cc        deltaz(k) = 2.d0*sin(hpin)*sin(dble(2*k-1)*hpin)
cc        deltaz(nz-k) = deltaz(k)
          deltaz(k)=abs(zets(k+1)-zets(k))
          deltaz(nz-k) = deltaz(k)
c          wzet(k) =abs(zets(k+1)-zets(k))
c          wzet(nz-k) = wzet(k)
ccccccccccccccccccccccccccccccccccccccccccccc GA
      enddo
      if (mod(nz,2).eq.0) deltaz(nz/2) = 2.d0*sin(hpin)*
     &  sin(dble(nzm)*hpin)
c
c Nadim: need to take care of below for tt2 field (svx and svy)
c     compute anisotropic spectral viscosity
c     x-direction
      do i=1,nx
         svx(i)=0.d0
      enddo
      if(epssvx.gt.0)then
         write(26,*)'Using spectral viscosity in x-direction'
         if(msvx.lt.nxh)then
            do i=msvx+1,nxh
               svx(i)=epssvx*exp((dble(i)**2-dble(nxh )**2)
     *                                   /(dble(i)**2-dble(msvx)**2))
               svx(nx+2-i)=svx(i)
            enddo
         else
          write(26,*)'Error in Spectral viscosity:  msvx >= nxh'
          write(*,*)'Error in Spectral viscosity:  msvx >= nxh'
          stop
         endif
      endif
c     y-direction
      do j=1,nyh
         svy(j)=0.d0
      enddo
      if(epssvy.gt.0)then
         write(26,*)'Using spectral viscosity in y-direction'
         if(msvy.lt.nyh)then
            do j=msvy+1,nyh
               svy(j)=epssvy*exp((dble(j)**2-dble(nyh )**2)
     *                                   /(dble(j)**2-dble(msvy)**2))
            enddo
         else
          write(26,*)'Error in Spectral viscosity:  msvy >= nyh'
          write(*,*)'Error in Spectral viscosity:  msvy >= nyh'
          stop
         endif
      endif
c
c     compute anisotropic spectral viscosity
c     x-direction
      do i=1,nx
         svx2(i)=0.d0
      enddo
      if(epssvx2.gt.0)then
         write(26,*)'Using spectral viscosity in x-direction for tt2'
         if(msvx2.lt.nxh)then
            do i=msvx2+1,nxh
               svx2(i)=epssvx2*exp((dble(i)**2-dble(nxh )**2)
     *                                   /(dble(i)**2-dble(msvx2)**2))
               svx2(nx+2-i)=svx2(i)
            enddo
         else
          write(26,*)'Error in Spectral viscosity:  msvx2 >= nxh'
          write(*,*)'Error in Spectral viscosity:  msvx2 >= nxh'
          stop
         endif
      endif
c     y-direction
      do j=1,nyh
         svy2(j)=0.d0
      enddo
      if(epssvy2.gt.0)then
         write(26,*)'Using spectral viscosity in y-direction for tt2'
         if(msvy2.lt.nyh)then
            do j=msvy2+1,nyh
               svy2(j)=epssvy*exp((dble(j)**2-dble(nyh )**2)
     *                                   /(dble(j)**2-dble(msvy2)**2))
            enddo
         else
          write(26,*)'Error in Spectral viscosity:  msvy2 >= nyh'
          write(*,*)'Error in Spectral viscosity:  msvy2 >= nyh'
          stop
         endif
      endif
c

c
c compute gravity vector 
cNOTE: 1) recall gy(nly0,nz) and gz(nyl0,nz) to de-aliase forcing
c      2) gz cannot be a function of y when substracting tt mean
c      3) recall that here zets are in [-1,1]
c
      FLGVARY=.false.
c
cNOTE: theta is now input in input file (chnl.ini)
c      theta=5.d0 !degrees
      theta=2.d0*pi/3.6d2*theta !rad
c
      if(FLGVARY)then
cNOTE: if FLGVARY=.true. check subroutine nltrms_up for 
c      settling velocity computation. Also check main.F
c      for particles velocities courant computation 
         do k=1,nz
            do j=1,nyl0
               gy(j,k)=0.d0 !dsin(theta)
               gz(j,k)=erf((1.d2*(zets(k)+1.d0)-1.d2)*3.d-1)!-dcos(theta)
            enddo
         enddo
      else !gy and gz do not vary with y -> only first value used
         do k=1,nz
c            gy(1,k)=0.d0                                  
c            gz(1,k)=erf((1.d2*(zets(k)+1.d0)-1.d2)*3.d-1)
c            gz(1,k)=0.d0
            gy(1,k)=dsin(theta) !1.d0
            gz(1,k)=-dcos(theta)!/dsin(theta)
         enddo
      endif
c
c---------------------------------------------------------------------
c       compute chb matrix for Chebyshev spectrum
c---------------------------------------------------------------------
c
c      pin=pi/dble(nzm) !pin is computed before for zets
      inzm=1.d0/dble(nzm)
      tinzm=2.d0*inzm
      hinzm=inzm/2.d0
C$OMP PARALLEL DEFAULT(SHARED)
C$OMP DO PRIVATE(i,j)
      do i=2,nzm
       do j=2,nzm
        chb(i,j)=tinzm*cos(pin*dble(j-1)*dble(i-1))
       enddo 
       chb(i,1) =inzm
       chb(i,nz)=inzm*cos(pi*dble(i-1))
      enddo
C$OMP END DO
C$OMP DO PRIVATE(j)
      do j=2,nzm
       chb(1,j) =inzm
       chb(nz,j)=inzm*cos(pi*dble(j-1))
      enddo
C$OMP END DO
C$OMP END PARALLEL
      chb(1,1)  =hinzm
      chb(1,nz) =hinzm
      chb(nz,1) =hinzm
      chb(nz,nz)=hinzm*cos(pi*dble(nzm))
c
c---------------------------------------------------------------------
c       compute sm_d1 matrix
c---------------------------------------------------------------------
c     ssc=0 will have spectral scheme in the z direction
      if(ssc.eq.0)then
        do j=2,nz
          do i=1,j-1
            chbd1(i,j) = ((-1)**(i+j))/
     &        (-2.d0*sin(hpin*dble(i-j))*sin(hpin*dble(i+j-2)))
          enddo
        enddo
c
c     fix edges
c
        call dscal(nzmm,2.d0, chbd1(1,2),nz0)
        call dscal(nzmm,0.5d0,chbd1(2,nz),1)
c
c     fill in the lower triangle
c
        do j=1,nzm
          do i=j+1,nz
            chbd1(i,j) = -chbd1(nz+1-i,nz+1-j)
          enddo
        enddo
c
c     compute the diagonal terms
c
        do i= 2,nz/2
          chbd1(i,i) = -0.5d0*zets(i)/sin(pin*dble(i-1))**2
          chbd1(nz+1-i,nz+1-i) = -chbd1(i,i)
        enddo
c
c     compute the (1,1) and (nz,nz) terms
c
        chbd1(1,1)   = (2.d0*dble(nzm**2)+1.d0)/6.d0
        chbd1(nz,nz) = -chbd1(1,1)
c
      write(100,*) 'z spectral'
      endif
ccccccccccccccccccccccccccccccccccccccccccccc GA ccccccccccccccccccccc
c--------------------------------------------------------
c     compute compact scheme d1 matrix
c--------------------------------------------------------
c     ssc=1 will have compact scheme in the z direction
      if(ssc.eq.1)then
      do j=1,nz
         do i=1,nz
            R1z(i,j) = 0.d0
            L1z(i,j) = 0.d0
         enddo
      enddo
c
      do k=1,nz
      if (k.eq.2) then
        aa=deltaz(k-1)
        bb=deltaz(k)
        dd=deltaz(k)+deltaz(k+1)
        ee=deltaz(k)+deltaz(k+1)+deltaz(k+2)
c
        a=(bb**2.d0*dd*ee*(6.d0*aa**3.d0+2.d0*bb*dd*ee+
     &     aa**2.d0*(4.d0*bb+5.d0*(dd+ee))+aa*(4.d0*dd*ee
     &     +3.d0*bb*(dd+ee))))/(aa*(aa+bb)**3.d0*(aa+dd)**2.d0
     &       *(aa+ee)**2.d0)
        b=(aa**2.d0*dd*ee*(aa*(-4.d0*bb**2.d0-2.d0*dd*ee
     &   +3.d0*bb*(dd+ee))+bb*(-6.d0*bb**2.d0-4.d0*dd*ee
     &    +5.d0*bb*(dd+ee))))
     &     /(bb*(aa+bb)**3.d0*(bb-dd)**2.d0*(bb-ee)**2.d0)
        c=-2.d0/aa+2.d0/bb+1.d0/dd+1.d0/ee
        d=(aa**2.d0*bb**2.d0*ee)/((bb-dd)**2.d0*dd
     &      *(aa+dd)**2.d0*(dd-ee))
        e=(aa**2.d0*bb**2.d0*dd)
     &     /((bb-ee)**2.d0*ee*(aa+ee)**2.d0*(-dd+ee))
        alfa=(bb**2.d0*dd*ee)/((aa+bb)**2.d0*(aa+dd)*(aa+ee))
        beta=(aa**2.d0*dd*ee)/((aa+bb)**2.d0*(bb-dd)*(bb-ee))
c        
        R1z(k,k-1)=a
        R1z(k,k)=c
        R1z(k,k+1)=b
        R1z(k,k+2)=d
        R1z(k,k+3)=e
c        
        L1z(k,k-1)=alfa
        L1z(k,k)=1.d0
        L1z(k,k+1)=beta
c            
      elseif (k.eq.(nz-1)) then
        aa=deltaz(k)
        bb=deltaz(k-1)
        dd=deltaz(k-1)+deltaz(k-2)
        ee=deltaz(k-1)+deltaz(k-2)+deltaz(k-3)
c        
        a=-(bb**2.d0*dd*ee*(6.d0*aa**3.d0+2.d0*bb*dd*ee
     &    +aa**2.d0*(4.d0*bb+5.d0*(dd+ee))+aa*(4.d0*dd*ee
     &     +3.d0*bb*(dd+ee))))/(aa*(aa+bb)**3.d0*
     &      (aa+dd)**2.d0*(aa+ee)**2.d0)
        b=-(aa**2.d0*dd*ee*(bb*(-6.d0*bb**2.d0-4.d0*dd*ee
     &    +5.d0*bb*(dd+ee))+aa*(-4.d0*bb**2.d0-2.d0*dd*ee
     &     +3.d0*bb*(dd+ee))))/(bb*(aa+bb)**3.d0*(bb-dd)**2.d0
     &      *(bb-ee)**2.d0)
        c=-(-2.d0/aa+2.d0/bb+1.d0/dd+1.d0/ee)
        d=-(aa**2.d0*bb**2.d0*ee)/((bb-dd)**2.d0*dd
     &     *(aa+dd)**2.d0*(dd-ee))
        e=-(aa**2.d0*bb**2.d0*dd)/((bb-ee)**2.d0*ee
     &     *(aa+ee)**2.d0*(-dd+ee))
        alfa=(aa**2.d0*dd*ee)/((aa+bb)**2.d0*(bb-dd)*(bb-ee))
        beta=(bb**2.d0*dd*ee)/((aa+bb)**2.d0*(aa+dd)*(aa+ee))
c
        R1z(k,k+1)=a
        R1z(k,k)=c
        R1z(k,k-1)=b
        R1z(k,k-2)=d
        R1z(k,k-3)=e
c        
        L1z(k,k-1)=alfa
        L1z(k,k)=1.d0
        L1z(k,k+1)=beta
c        
      elseif (k.eq.1) then
        bb=deltaz(k)
        cc=deltaz(k)+deltaz(k+1)
        dd=deltaz(k)+deltaz(k+1)+deltaz(k+2)
        ee=deltaz(k)+deltaz(k+1)+deltaz(k+2)+deltaz(k+3)
        ff=deltaz(k)+deltaz(k+1)+deltaz(k+2)+deltaz(k+3)
     &         +deltaz(k+4)
c        
        alfa=cc*dd*ee*ff/((cc-bb)*(dd-bb)*(ee-bb)*(ff-bb))
        a=1.d0/cc+1.d0/dd+1.d0/ee+1.d0/ff+2.d0/bb
        b=alfa*(1.d0/(cc-bb)+1.d0/(dd-bb)+1.d0/(ee-bb)
     &                              +1.d0/(ff-bb)-2.d0/bb)
        c=(bb/(cc-bb))**2.d0*dd*ee*ff/((dd-cc)*(ee-cc)*(ff-cc)*(-cc))
        d=(bb/(bb-dd))**2.d0*cc*ee*ff/((cc-dd)*(ee-dd)*(ff-dd)*(-dd))
        e=(bb/(bb-ee))**2.d0*cc*dd*ff/((cc-ee)*(dd-ee)*(ff-ee)*(-ee))
        f=(bb/(bb-ff))**2.d0*cc*dd*ee/((cc-ff)*(dd-ff)*(ee-ff)*(-ff))
c        
        R1z(k,k)=a
        R1z(k,k+1)=b
        R1z(k,k+2)=c 
        R1z(k,k+3)=d
        R1z(k,k+4)=e
        R1z(k,k+5)=f
c        
        L1z(k,k)=1.d0
        L1z(k,k+1)=alfa
c        
      elseif (k.eq.nz) then
        bb=deltaz(k-1)
        cc=deltaz(k-1)+deltaz(k-2)
        dd=deltaz(k-1)+deltaz(k-2)+deltaz(k-3)
        ee=deltaz(k-1)+deltaz(k-2)+deltaz(k-3)+deltaz(k-4)
        ff=deltaz(k-1)+deltaz(k-2)+deltaz(k-3)+deltaz(k-4)
     &                         +deltaz(k-5)
c        
        alfa=cc*dd*ee*ff/((cc-bb)*(dd-bb)*(ee-bb)*(ff-bb))
        a=-1.d0/cc-1.d0/dd-1.d0/ee-1.d0/ff-2.d0/bb
        b=alfa*(-1.d0/(cc-bb)-1.d0/(dd-bb)-1.d0/(ee-bb)
     &                             -1.d0/(ff-bb)+2.d0/bb)
        c=(-bb/(bb-cc))**2.d0*(-dd)*(-ee)*(-ff)/((cc-dd)
     &                             *(cc-ee)*(cc-ff)*(cc))
        d=(-bb/(dd-bb))**2.d0*(-cc)*(-ee)*(-ff)/((dd-cc)
     &                             *(dd-ee)*(dd-ff)*(dd))
        e=(-bb/(ee-bb))**2.d0*(-cc)*(-dd)*(-ff)/((ee-cc)
     &                             *(ee-dd)*(ee-ff)*(ee))
        f=(-bb/(ff-bb))**2.d0*(-cc)*(-dd)*(-ee)/((ff-cc)
     &                             *(ff-dd)*(ff-ee)*(ff))
c        
        R1z(k,k)=a
        R1z(k,k-1)=b
        R1z(k,k-2)=c 
        R1z(k,k-3)=d
        R1z(k,k-4)=e
        R1z(k,k-5)=f
c        
        L1z(k,k)=1.d0
        L1z(k,k-1)=alfa
c        
      else 
        aa=deltaz(k-1)+deltaz(k-2)
        bb=deltaz(k-1)
        dd=deltaz(k)
        ee=deltaz(k)+deltaz(k+1)
c
        a=(bb**2.d0*dd**2.d0*ee)/(aa*(aa-bb)**2.d0
     &     *(aa+dd)**2.d0*(aa+ee))
        b=(aa*dd**2.d0*ee*(aa*(5.d0*bb**2.d0+3.d0*bb*dd
     &    +4.d0*bb*ee+2.d0*dd*ee)-bb*(6.d0*bb**2.d0+4.d0*bb*dd
     &     +5.d0*bb*ee+3.d0*dd*ee)))/((aa-bb)**2.d0*bb*
     &      (bb+dd)**3.d0*(bb+ee)**2.d0)
        c=-1.d0/aa-2.d0/bb+2.d0/dd+1.d0/ee
        d=(aa*bb**2.d0*ee*(aa*dd*(3.d0*bb+5.d0*dd)
     &    -2.d0*aa*(bb+2.d0*dd)*ee+dd*(4.d0*bb*dd
     &     +6.d0*dd**2.d0-3.d0*bb*ee-5.d0*dd*ee)))
     &      /(dd*(aa+dd)**2.d0*(bb+dd)**3.d0*(dd-ee)**2.d0)
        e=-(aa*bb**2.d0*dd**2.d0)/((dd-ee)**2.d0*ee*(aa+ee)
     &     *(bb+ee)**2.d0)
        alfa=(aa*dd**2.d0*ee)/((aa-bb)*(bb+dd)**2.d0*(bb+ee))
        beta=(aa*bb**2.d0*ee)/((aa+dd)*(bb+dd)**2.d0*(-dd+ee))
c        
          R1z(k,k-2)=a
          R1z(k,k-1)=b
          R1z(k,k)=c
          R1z(k,k+1)=d
          R1z(k,k+2)=e
c        
          L1z(k,k-1)=alfa
          L1z(k,k)=1.d0
          L1z(k,k+1)=beta
        endif
       enddo
c
      call dgetrf(nz,nz,L1z,nz,ipiv1,i22)
      call dgetri(nz,L1z,nz,ipiv1,work,lwork,i33)
c
      call dgemm('N','N',nz,nz,nz,1.d0,L1z,nz,
     &                   R1z,nz,0.d0,chbd1,nz0)
c
      write(100,*) 'z compact'
      endif
      close(unit=100)
ccccccccccccccccccccccccccccccccccccccccccccc GA
c---------------------------------------------------------------------
c       stretch for Lz.neq.2: z->(-rlenz/2,rlenz/2)
c---------------------------------------------------------------------
      if(rlenz.ne.2.d0)then
c     stretch grid. nothing to be done to weights
         do k=1,nz
            zets(k)  = rlenz/2.d0*(zets(k)-1)+rlenz/2.d0
c            wzet(k)  = wzet(k)*rlenz !*rlenz/2.d0 !2.d0/rlenz
            if(k.lt.nz) deltaz(k)= deltaz(k)*rlenz/2.d0
         enddo
c     
c     stretch chbd1. NOTE: chbd2 automatically done
         call dscal(nz0z,2.d0/rlenz,chbd1,1)
      endif
c
      do k=1,nzm
         ndeltaz(k) = -deltaz(k)
         idz(k) = 1.d0/deltaz(k)
         nidz(k) = -idz(k)
      enddo
c     
      write(26,119)
      do k = 1,nzm
         write(26,121) k,zets(k),wzet(k),deltaz(k)
      enddo
      write(26,121) nz,zets(nz),wzet(k)
 119  FORMAT(/,'   i   gauss lobatto pts     weights         spacings')
 121  FORMAT(2x,i4,3(2x,1pe15.7,2x,1pe15.7,2x,1pe15.7))
c
c     compute locations to output time series
      k=nz
      do j=1,11
        if(kts(j).ge.re) then
          kts(j)=1
        else
          k=k-1
          do while((rlenz/2.d0+zets(k))*re-dble(kts(j))
     &            .lt.0.d0)
             k=k-1
          enddo
          kts(j)=k
        endif
c       write(*,*) j, kts(j)
      enddo
c 
ccccccccccccccccccccccccccccccccccccccccccccccc GA ccccccccccccccccccc    
c---------------------------------------------------------------------
c     Parameters for the Immersed Boundary Method
c---------------------------------------------------------------------
      if (FLIBM) then 
c
         do k=1,nz
            do j=1,nz
               chbd1_ibm(j,k)=chbd1(j,k)
            enddo
         enddo
         do k=1,nz
            chbd1_ibm(1,k)=0
         enddo
         chbd1_ibm(1,1)=1
         call dgetrf(nz,nz,chbd1_ibm,nz,ipiv1,i22)
         call dgetri(nz,chbd1_ibm,nz,ipiv1,work,lwork,i33)
c
         do n=1,num_p
            do l=1,n_ll
               x_o        (l,n)=0.d0
               y_o        (l,n)=0.d0
               z_o        (l,n)=0.d0
               dv_l       (l,n)=0.d0
               u_p        (l,n)=0.d0
               v_p        (l,n)=0.d0
               w_p        (l,n)=0.d0
            enddo
         enddo
         z_o_avg = 0
c
c Nadim there should be an if(FLTHRM) for tt_p, tt2_p, ...  above
      endif
ccccccccccccccccccccccccccccccccccccccccccccccc GA
c---------------------------------------------------------------------
c     compute sm_d2 matrix by forming the product sm_d1*sm_d1
c---------------------------------------------------------------------
c       
      call dgemm('N','N',nz,nz,nz,1.d0,chbd1,nz0,chbd1,nz0,0.d0,
     &  chbd2,nz0)
c
c---------------------------------------------------------------------
c     compute sm_d2 matrix with dirichlet-dirichlet bc's
c---------------------------------------------------------------------
c
      do j = 1,nz
        do i = 2,nzm
          ci_drdr(i,j) = chbd2(i,j)
        enddo
      enddo
      ci_drdr(1,1)   = 1.d0
      ci_drdr(nz,nz) = 1.d0
c
c---------------------------------------------------------------------
c     diagonalize the ci_drdr matrix 
c---------------------------------------------------------------------
c
      write(26,*)' 1- & infinity-norm condition numbers, ',
     &  'max and min eigenvalues'
      call e_decomp(nz,ci_drdr,nz0,cn_drdr,e_drdr,condo,condi,emax,emin)
      write(26,*)'cn_drdr: ', condo, condi, emax, emin
c
c---------------------------------------------------------------------
c     compute the sm_d2 matrix with neumann-dirichlet bc's
c---------------------------------------------------------------------
c
      do j = 1,nz
        ci_nodr(1,j)  = chbd1(1,j)
        do i = 2,nzm
          ci_nodr(i,j) = chbd2(i,j)
        enddo
      enddo
      ci_nodr(nz,nz) = 1.d0
c
c---------------------------------------------------------------------
c     diagonalize the ci_nodr matrix
c---------------------------------------------------------------------
c
      call e_decomp(nz,ci_nodr,nz0,cn_nodr,e_nodr,condo,condi,emax,emin)
      write(26,*)'cn_nodr: ', condo, condi, emax, emin
c
c---------------------------------------------------------------------
c     compute the v_d2 matrix with mixed-mixed bc's
c---------------------------------------------------------------------
c uvbb and uvbt are input in input file
      uvat=1.d0-uvbt
      uvab=1.d0-uvbb
c in case we need to add some shear at the walls modify these accordingly
      ugt=0.d0 !homogeneous bc at top    wall
      ugb=0.d0 !homogeneous bc at bottom wall
      vgt=0.d0 !homogeneous bc at top    wall
      vgb=0.d0 !homogeneous bc at bottom wall
c
      d00   = uvat + uvbt *chbd1(1,1)
      d0n   = uvbt * chbd1(1,nz)
      dn0   = uvbb * chbd1(nz,1)
      dnn   = uvab + uvbb * chbd1(nz,nz)
      iden  = 1.d0/(dn0*d0n-d00*dnn)
      do i=1,nzmm
        uvprbc1(i) = (dnn*(uvbt * chbd1( 1,i+1))
     &              -d0n*(uvbb * chbd1(nz,i+1))) * iden
        uvprbcn(i) = (d00*(uvbb * chbd1(nz,i+1))
     &              -dn0*(uvbt * chbd1(1, i+1))) * iden
      enddo
c
      uvmbct1=dnn*(-iden)
      uvmbcb1=d0n*(-iden)
      uvmbctn=dn0*(-iden)
      uvmbcbn=d00*(-iden)
c
      do j = 2,nzm
        do i = 2,nzm
          uvci_mimi(i-1,j-1) = chbd2(i,j) + chbd2(i, 1)*uvprbc1(j-1)
     &                                    + chbd2(i,nz)*uvprbcn(j-1)
        enddo
      enddo
c---------------------------------------------------------------------
c     diagonalize the vci_mimi matrix
c---------------------------------------------------------------------
c
      call e_decomp(nzmm,uvci_mimi,nzmm0,uvcn_mimi,uve_mimi,condo,condi,
     &  emax,emin)
      write(26,*)'uvcn_mimi: ', condo, condi, emax, emin
c
c---------------------------------------------------------------------
c     compute the w_d2 matrix with mixed-mixed bc's
c---------------------------------------------------------------------
c these are values for homogeneous dirichlet-dirichlet bc's 
      wat=1.d0
      wab=1.d0
      wbt=0.d0
      wbb=0.d0
      wgt=0.d0 !no penetration at top    wall
      wgb=0.d0 !no penetration at bottom wall
c
      d00   = wat + wbt *chbd1(1,1)
      d0n   = wbt * chbd1(1,nz)
      dn0   = wbb * chbd1(nz,1)
      dnn   = wab + wbb * chbd1(nz,nz)
      iden  = 1.d0/(dn0*d0n-d00*dnn)
      do i=1,nzmm
        wprbc1(i) = (dnn*(wbt * chbd1( 1,i+1))
     &              -d0n*(wbb * chbd1(nz,i+1))) * iden
        wprbcn(i) = (d00*(wbb * chbd1(nz,i+1))
     &              -dn0*(wbt * chbd1(1, i+1))) * iden
      enddo
c
      wmbct1=dnn*(-iden)
      wmbcb1=d0n*(-iden)
      wmbctn=dn0*(-iden)
      wmbcbn=d00*(-iden)
c
      do j = 2,nzm
        do i = 2,nzm
          wci_mimi(i-1,j-1) = chbd2(i,j) + chbd2(i, 1)*wprbc1(j-1)
     &                                   + chbd2(i,nz)*wprbcn(j-1)
        enddo
      enddo
c---------------------------------------------------------------------
c     diagonalize the wci_mimi matrix
c---------------------------------------------------------------------
      call e_decomp(nzmm,wci_mimi,nzmm0,wcn_mimi,we_mimi,condo,condi,
     &  emax,emin)
      write(26,*)'wcn_mimi: ', condo, condi, emax, emin
c
c---------------------------------------------------------------------
c     compute the tt_d2 matrix with mixed-mixed bc's
c---------------------------------------------------------------------
c
      if(FLTHRM)then
cccc
       if(ttat(1).eq.999)then
        ttat(1)=-re*pr*ws*gz(1,1)/sqrt(gy(1,1)**2+gz(1,1)**2)
        write(26,*)'ttat(1) computed internally'
        write( *,*)'ttat(1) computed internally'
       endif
       if(ttab(1).eq.999)then
        ttab(1)=-re*pr*ws*gz(1,nz)/sqrt(gy(1,nz)**2+gz(1,nz)**2)
        write(26,*)'ttab(1) computed internally'
        write( *,*)'ttat(1) computed internally'
       endif
cccc
      d00   = ttat(1) + ttbt(1) *chbd1(1,1)
      d0n   = ttbt(1) * chbd1(1,nz)
      dn0   = ttbb(1) * chbd1(nz,1)
      dnn   = ttab(1) + ttbb(1) * chbd1(nz,nz)
      iden  = 1.d0/(dn0*d0n-d00*dnn)
      do i=1,nzmm
        ttprbc1(i) = (dnn*(ttbt(1) * chbd1( 1,i+1))
     &              -d0n*(ttbb(1) * chbd1(nz,i+1))) * iden
        ttprbcn(i) = (d00*(ttbb(1) * chbd1(nz,i+1))
     &              -dn0*(ttbt(1) * chbd1(1, i+1))) * iden
      enddo
c
      ttmbct1=dnn*(-iden)
      ttmbcb1=d0n*(-iden)
      ttmbctn=dn0*(-iden)
      ttmbcbn=d00*(-iden)
c
      do j = 2,nzm
        do i = 2,nzm
          ttci_mimi(i-1,j-1) = chbd2(i,j) + chbd2(i, 1)*ttprbc1(j-1)
     &                                    + chbd2(i,nz)*ttprbcn(j-1)
        enddo
      enddo
c---------------------------------------------------------------------
c     diagonalize the ttci_mimi matrix
c---------------------------------------------------------------------
      call e_decomp(nzmm,ttci_mimi,nzmm0,ttcn_mimi,tte_mimi,condo,condi,
     &  emax,emin)
      write(26,*)'ttcn_mimi: ', condo, condi, emax, emin
c
      endif
c
      if(FLTHRM2)then
cccc
       if(tt2at(1).eq.999)then
        tt2at(1)=-re*pr2*ws*gz(1,1)/sqrt(gy(1,1)**2+gz(1,1)**2)
        write(26,*)'tt2at(1) computed internally'
        write( *,*)'tt2at(1) computed internally'
       endif
       if(tt2ab(1).eq.999)then
        tt2ab(1)=-re*pr2*ws*gz(1,nz)/sqrt(gy(1,nz)**2+gz(1,nz)**2)
        write(26,*)'tt2ab(1) computed internally'
        write( *,*)'tt2at(1) computed internally'
       endif
cccc
      d00   = tt2at(1) + tt2bt(1) *chbd1(1,1)
      d0n   = tt2bt(1) * chbd1(1,nz)
      dn0   = tt2bb(1) * chbd1(nz,1)
      dnn   = tt2ab(1) + tt2bb(1) * chbd1(nz,nz)
      iden  = 1.d0/(dn0*d0n-d00*dnn)
      do i=1,nzmm
        tt2prbc1(i) = (dnn*(tt2bt(1) * chbd1( 1,i+1))
     &              -d0n*(tt2bb(1)   * chbd1(nz,i+1))) * iden
        tt2prbcn(i) = (d00*(tt2bb(1) * chbd1(nz,i+1))
     &              -dn0*(tt2bt(1)   * chbd1(1, i+1))) * iden
      enddo
c
      tt2mbct1=dnn*(-iden)
      tt2mbcb1=d0n*(-iden)
      tt2mbctn=dn0*(-iden)
      tt2mbcbn=d00*(-iden)
c
      do j = 2,nzm
        do i = 2,nzm
          tt2ci_mimi(i-1,j-1) = chbd2(i,j) + chbd2(i, 1)*tt2prbc1(j-1)
     &                                     + chbd2(i,nz)*tt2prbcn(j-1)
        enddo
      enddo
c---------------------------------------------------------------------
c     diagonalize the ttci_mimi matrix
c---------------------------------------------------------------------
      call e_decomp(nzmm,tt2ci_mimi,nzmm0,tt2cn_mimi,tt2e_mimi,condo,
     &  condi,emax,emin)
      write(26,*)'tt2cn_mimi: ', condo, condi, emax, emin
c
      endif
c---------------------------------------------------------------------
c     compute the sm_d2 matrix with neumann-neumann bc's 
c---------------------------------------------------------------------
c
      d00   = chbd1(1,1)
      d0n   = chbd1(1,nz)
      dn0   = chbd1(nz,1)
      dnn   = chbd1(nz,nz)
      iden  = 1.d0/(dn0*d0n-d00*dnn)
      do i=1,nzmm
        prbc1(i) = (dnn*chbd1(1, i+1)-d0n*chbd1(nz,i+1))*iden
        prbcn(i) = (d00*chbd1(nz,i+1)-dn0*chbd1(1, i+1))*iden
      enddo
c
      do j = 2,nzm
        do i = 2,nzm
          ci_nono(i-1,j-1) = chbd2(i,j) + chbd2(i, 1)*prbc1(j-1)
     &                                  + chbd2(i,nz)*prbcn(j-1)
        enddo
      enddo
c
c---------------------------------------------------------------------
c     diagonalize the ci_nono matrix
c---------------------------------------------------------------------
c
      call e_decomp(nzmm,ci_nono,nzmm0,cn_nono,e_nono,condo,condi,
     &  emax,emin)
      write(26,*)'cn_nono: ', condo, condi, emax, emin
c
c---------------------------------------------------------------------
c     read in data
c---------------------------------------------------------------------
c
      rtime = 0.d0
      if (ch_fin.eq.0) then
c     initialize field for zero start run for ch_fin = 0
c        call init_uvw
c        call courant(dt0,u,v,w)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do i=1,nx0
         do j=1,ny0
            do k=1,nz
               u(i,j,k)=0.d0
c               v(i,j,k)= 0.5d0*re*(-ppA)*(zets(k))**2 
c     &           -0.5d0*re*(-ppA)*rlenz*(zets(k))
c     &           -0.125d0*(-ppA)*re*(2.d0*(amp+gap)-3.d0*rlenz+
c     &              2.d0*amp*dcos(2.d0*pi*dble(j-1)*deltay/period))*
c     &             (2.d0*(amp+gap)-rlenz+ 
c     &              2.d0*amp*dcos(2.d0*pi*dble(j-1)*deltay/period))
               v(i,j,k)=1.d0
               w(i,j,k)=0.d0
               tt(i,j,k)=1.d0
               tt2(i,j,k)=1.0d0
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
c
c      thickness = 1.d-2
c      center_tt = rleny/2.d0
c      hwidth_tt = 8.5d0
ccC$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
c      do k=1,nz
c         do j=1,ny0
c            do i=1,nx0
cc               if(zets(k) .lt. (zets(nz)+0.2)) then
cc                  goto 114
c               if(yets(j) .ge. center_tt) then
c                  tt(i,j,k)=0.5d0-0.5d0*dtanh((yets(j)-
c     &                             (center_tt+hwidth_tt))/thickness)
c               else
c                  tt(i,j,k)=0.5d0+0.5d0*dtanh((yets(j)-
c     &                             (center_tt-hwidth_tt))/thickness)
c               endif
cc               endif
c            enddo
c         enddo
c      enddo
c      do k=1,nz
c         do j=1,ny0
c            do i=1,nx0
c               if(zets(k) .gt. 0) then
c                  tt(i,j,k)=0.d0
c               endif
c            enddo
c         enddo
c      enddo
cC$OMP END PARALLEL DO
c      do k=k_ibm+1,nz
c         do j=1,ny0
c            do i=1,nx0
c               tt(i,j,k)=0.d0
c            enddo
c         enddo
c      enddo
c
         call initial_lag
         call lagrangian_pts_gen
         if(itime.eq.0) call output_2d
      else if (ch_fin.eq.1) then
c     dt_old is read from input file
c     read field data as written by routine "output" for ch_fin = 1
        call initial_lag
        call input
        call lagrangian_pts_gen
        if(itime.eq.0) call output_2d
      else
        write(26,*) 'invalid value for ch_fin'
        stop
      endif
c
      if(m_fin.eq.0)then
         do i=1,np
            NT(i)=0
         enddo
      elseif(m_fin.eq.1)then
c     read mean value and NT
         call input_means
      else
         write(26,*) 'invalid value for m_fin'
         stop
      endif
c
      if (FLTIMER) call end_timer(10)
      return
122     END
c
c
c
c     Eigendecomposition
c
c     Input:  the matrix a
c     Output: eval contains the eigenvalues of a
c             evec contains the eigenvector matrix of a
c             a    now contains the inverse of evec
c
c     It is assumed that n <= nz and n0 <= nz0
c
      SUBROUTINE e_decomp(n,a,n0,evec,eval,condo,condi,evalmax,evalmin)
#include "common.inc"
      INTEGER n,n0
      REAL*8 a(n0,n),evec(n0,n),eval(n),condo,condi,evalmax,evalmin
c
      INTEGER nwork
      REAL*8 eps
      PARAMETER (nwork=nz*(nz+4),eps=1.d-9)
c
      INTEGER i0,i1,i2,i3,i,j,ipiv(nz)
      REAL*8 val,evi(nz),work(nwork)
#ifdef ESSL
c for ESSL
      COMPLEX*16 ceval(nz),cevec(nz0,nz)
      LOGICAL alltrue(nz)
      REAL*8 dnorm2
      EXTERNAL dnorm2
#else
      REAL*8 dnrm2
      EXTERNAL dnrm2
#endif
c      INTEGER idamax,idmax,idmin
      INTEGER idamax,auxidmax,auxidmin
      REAL*8 dasum
c      EXTERNAL idamax,idmax,idmin,dasum
      EXTERNAL idamax,auxidmin,auxidmax,dasum
c
      i0=0
#ifdef ESSL
c for ESSL
      do i=1,n
        alltrue(i) = .true.
      enddo
      call dgeev(2,a,n0,ceval,cevec,nz0,alltrue,n,work,nwork)
      do i=1,n
        eval(i) = dreal(ceval(i))
        evi(i)  = dimag(ceval(i))
      enddo
      do j=1,n
        do i=1,n
          evec(i,j) = dreal(cevec(i,j))
        enddo
        val = 1.d0/dnorm2(n,evec(1,j),1)
        call dscal(n,val,evec(1,j),1)
      enddo
#else
      call dgeev('N','V',n,a,n0,eval,evi,evec,n0,evec,n0,work,nwork,i0)
      do j=1,n
        val = 1.d0/dnrm2(n,evec(1,j),1)
        call dscal(n,val,evec(1,j),1)
      enddo
#endif
c
      i1 = 0
      val = abs(evi(idamax(n,evi,1)))
      if (val.gt.eps) i1 = 1
      write(26,*) 'val, eps:',val,eps
c
c      evalmax = eval(idmax(n,eval,1))
      evalmax = eval(auxidmax(n,eval,1))
c      evalmin = eval(idmin(n,eval,1))
      evalmin = eval(auxidmin(n,eval,1))
c
      do i=1,n
        work(i) = dasum(n,evec(1,i), 1)
        evi( i) = dasum(n,evec(i,1),n0)
      enddo
c      condo = work(idmax(n,work,1))
      condo = work(auxidmax(n,work,1))
c      condi = evi( idmax(n, evi,1))
      condi = evi( auxidmax(n, evi,1))
c
      call dcopy(n*n0,evec,1,a,1)
      call dgetrf(n,n,a,n0,ipiv,i2)
      call dgetri(n,a,n0,ipiv,work,nwork,i3)
c
      do i=1,n
        work(i) = dasum(n,a(1,i), 1)
        evi( i) = dasum(n,a(i,1),n0)
      enddo
c      condo = condo*work(idmax(n,work,1))
      condo = condo*work(auxidmax(n,work,1))
c      condi = condi*evi( idmax(n, evi,1))
      condi = condi*evi( auxidmax(n, evi,1))
c
      if (i0.eq.0.and.i1.eq.0.and.i2.eq.0.and.i3.eq.0) return
      write(26,*) 'error in e_decomp'
      write(26,*) 'i0, i1, i2, i3:',i0,i1,i2,i3
      print *, 'error in e_decomp'
c      stop
      END
c
c
c
c     Bad for parallel
c
      SUBROUTINE init_uvw
#include "common.inc"
c
      INTEGER i,j,k
      INTEGER jran
      REAL*8  rnum
      REAL*8  fac,sum,sum1
      INTEGER idum,i2d
      REAL*8  dum
cccccccccccccccccccccccccccccccccccccccccc GA cccccccccccc
      INTEGER nx0o,ny0o,nzo,itimeo !ICvel.ini
      REAL*8  rtimeo,dto,rlenzo,deltaxo,deltayo !ICvel.ini 
      REAL*8  ut(96,96,97),vt(96,96,97),wt(96,96,97)
      INTEGER nzdft,nzdftm,l,flag
      REAL*8 pindft,hpindft,deltazdft(97),zetsdft(97),r,rr(nz)
cccccccccccccccccccccccccccccccccccccccccc GA
c
      fac=1.d-2
      i2d=0
cccccccccccccccccccccccccccccccccccccccccc
c This is to input from an initialization file
c set .false. to .true.
      if(.true.)then
c         open(11,file='ICvel.ini',form=format_mode)
         open(11,file='ICvel.ini',form=format_mode,
     &       access='direct',
     &       recl=22500208)
          read(11,rec=1) idum,rtimeo,itimeo,dto,nx0o,ny0o,
     &                  nzo,deltaxo,deltayo,rlenzo,idum,idum,u,v,w,
     &                  idum,idum,pbxnt,pcxnt,pbynt,pcynt,pbxot,
     &                  pcxot,pbyot,pcyot,idum,idum,
     &                  pbxnb,pcxnb,pbynb,pcynb,pbxob,
     &                  pcxob,pbyob,pcyob,idum
         close(11)
         write(*,*)'nx0o:',nx0o,'  ny0o:',nx0o,'  nzo:',nzo
         write(*,*)'dto:',dto,'  rtimeo:',rtimeo,'  itimeo:',itimeo
         write(*,*)'deltaxo:',deltaxo,'deltyo:',deltayo,'rlenzo:',rlenzo
         write(*,*)'maxval(u):',maxval(u),'maxval(v):',maxval(v),
     $             'maxval(w):',maxval(w)
c         stop

c         print *,format_mode
c         rewind(11)
c         print *,'before'
cc         read(11) dum,idum,dum,idum,idum,idum,dum,dum,dum
c         read(11) dum
c         print *,dum
c         read(11) idum
c         print *,idum
c         read(11) dum
c         print *,dum
c         read(11) idum
c         print *,idum
c         read(11) idum
c         print *,idum
c         read(11) idum
c         print *,idum
c         read(11) dum
c         print *,dum
c         read(11) dum
c         print *,dum
c         read(11) dum
c         print *,dum
c         read(11) u,v,w
c         print *,'after2'
cc         read(11) pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
cc         read(11) pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
c         print *,'end'
c         close(11)
cccccccccccccccccccccccccccccccccccccccccc GA cccccccccccc
      goto 117
      nzdft = 97
      nzdftm = nzdft-1
      pindft  = pi/dble(nzdftm)
      hpindft = 0.5d0*pindft
      do j = 1,nzdft/2
        zetsdft(j) = cos(dble(j-1)*pindft)
        zetsdft(nzdft+1-j) = -zetsdft(j)
      enddo
c
      do k = 1,nzdftm/2
        deltazdft(k) = 2.d0*sin(hpindft)*sin(dble(2*k-1)*hpindft)
        deltazdft(nz-k) = deltazdft(k)
      enddo
c
      do k=1,nzdft
         zetsdft(k)  = rlenz/2.d0*(zetsdft(k)-1)+rlenz/2.d0
         deltazdft(k)= deltazdft(k)*rlenz/2.d0
      enddo 
c
      do i=1,96
         do j=1,96
            u(i,j,1)=ut(i,j,1)
            u(i,j,nz)=ut(i,j,nzdft)
         enddo
      enddo
c
      do k=1,nz
        rr(k)=0
      enddo
      do i=1,96
         do j=1,96
            do k=1,nz
               flag=1
               do l=1,nzdft
                  if(zets(k).le.zetsdft(nzdft+1-l).and.flag.eq.1)then
                     r=abs((zets(k)-zetsdft(nzdft+1-l))/
     &                   (zetsdft(nzdft+1-l+1)-zetsdft(nzdft+1-l)))
                     if(i.eq.45.and.j.eq.45)then
                        rr(k)=r
                     endif
                     u(i,j,k)=r*(ut(i,j,nzdft+1-l+1)-ut(i,j,nzdft+1-l))
     &                          +ut(i,j,nzdft+1-l)
                     v(i,j,k)=r*(vt(i,j,nzdft+1-l+1)-vt(i,j,nzdft+1-l))
     &                          +vt(i,j,nzdft+1-l)
                     w(i,j,k)=r*(wt(i,j,nzdft+1-l+1)-wt(i,j,nzdft+1-l))
     &                          +wt(i,j,nzdft+1-l)
                     flag=-1
                  endif
               enddo
            enddo
         enddo
      enddo
117   continue
cccccccccccccccccccccccccccccccccccccccc GA
c
         do k=1,nz
            do i=1,nx0
               u(i,ny0,k)=0
               v(i,ny0,k)=0
               w(i,ny0,k)=0
            enddo
            do j=1,ny
               u(nx0,j,k)=0
               v(nx0,j,k)=0
               w(nx0,j,k)=0
            enddo
         enddo
c
         if(FLTHRM)then
            open(11,file='ICconc.ini',form=format_mode)
            rewind(11)
            read(11) dum,idum,dum,idum,idum,idum,dum,dum,dum
            read(11) tt
            close(11)
c     This is no normalize tt to integrate to 1
            sum =0.d0
            sum1 =0.d0
            do k=1,nz
               do j=1,ny
                  do i=1,nx
                     call randomgen(jran,rnum)
                     tt(i,j,k)=(1.d0+fac*rnum)
                     sum=sum+tt(i,j,k)*wzet(k)
                     sum1=sum1+wzet(k)
                  enddo
               enddo
            enddo
            sum=sum/dble(nxny)
            write(*,*)'sum:',sum
            write(*,*)'sum1:',sum1/dble(nxny)
            sum1=0.d0
            do k=1,nz
               do j=1,ny
                  do i=1,nx
                     tt(i,j,k)=tt(i,j,k)/sum
                     sum1=sum1+tt(i,j,k)*wzet(k)
                  enddo
               enddo
            enddo
            write(*,*)'int(tt)/(Lx*Ly*Lz):',sum1/dble(nxny)
         endif
      endif
cccccccccccccccccccccccccccccccccccccccccc
cc NOTE u, v, and w are initialized to zero above, no need to repeat
      if(.false..and.FLTHRM)then
         sum =0.d0
         do k=1,nz
            do j=1,ny
c               if(i2d.eq.1)call randomgen(jran,rnum)
               do i=1,nx
c                  if(i2d.eq.0)call randomgen(jran,rnum)
c                  tt(i,j,k)=(1.d0+fac*rnum)
                  sum=sum+tt(i,j,k)*wzet(k)
               enddo
            enddo
         enddo         
         sum=sum/dble(nxny)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  tt(i,j,k)=tt(i,j,k)/sum
               enddo
            enddo
         enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c
      call fft_r2f_2d_new(u) 
      call fft_r2f_2d_new(v) 
      call fft_r2f_2d_new(w) 
      if(FLTHRM )call fft_r2f_2d_new(tt)
      if(FLTHRM2)call fft_r2f_2d_new(tt2)
      do k=1,nz
         u(1,1,k)=0.d0
         w(1,1,k)=0.d0
      enddo
      write(16,111)(zets(k),v(1,1,k),u(1,1,k),w(1,1,k),tt(1,1,k),
     &       tt2(1,1,k),k=1,nz)
      call fft_f2r_2d_new(u) 
      call fft_f2r_2d_new(v) 
      call fft_f2r_2d_new(w) 
      if(FLTHRM )call fft_f2r_2d_new(tt)
      if(FLTHRM2)call fft_f2r_2d_new(tt2)
 111  FORMAT(f12.8,5x,f12.8,5x,f12.8,5x,f12.8,5x,f12.8)
c
      return
      END
c
c This randomgen subroutine can be changed with any other subroutine we want
      subroutine randomgen(iseed,rnum)
      integer iseed
      real*8 rnum
      iseed  = mod(iseed*1366+150889,714025)
      rnum   = 2.d0*(dble(iseed)/714025.d0) - 1.d0
      return
      end

