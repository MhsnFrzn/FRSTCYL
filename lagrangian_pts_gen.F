      SUBROUTINE initial_lag
#include "common.inc"
c
      INTEGER i,j,k,n,l
      REAL*8 x2d(nx0,ny0), y2d(nx0,ny0),vertical
c Variables for scour Nadim
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      REAL*8 x01, y01, x02, y02, fac
      REAL*8 dumx, dumy, fac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c      call input_2d
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          x2d(i,j) = xets(i)
          y2d(i,j) = yets(j)
c          z2d(i,j) = rlenz/2.d0 -1.d0
        enddo
      enddo
C$OMP END PARALLEL DO
c
c      h0_sum = (rlenz/2.d0-1.d0)*dble(nx*ny)

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          x_o((i-1)*ny+j,1)  = xets(i)
          y_o((i-1)*ny+j,1)  = yets(j)
          z_o((i-1)*ny+j,1)  =-1.d0

          x_o((i-1)*ny+j,2)  = xets(i)
          y_o((i-1)*ny+j,2)  = yets(j)

          x_o((i-1)*ny+j,4)  = xets(i)
          y_o((i-1)*ny+j,4)  = yets(j)

          x_o((i-1)*ny+j,3)  = xets(i)
          y_o((i-1)*ny+j,3)  = yets(j)
          z_o((i-1)*ny+j,3)  = 1.d0
        enddo
      enddo
C$OMP END PARALLEL DO
c
cc Nadim Scour
      Rpillar = 2.5d0
c      nzp = nz-6
      fac = 0.8 ! spacing along circumference is less than deltay
      ! This is to account for radius of tube increasing over time
c      nxyp = int(ceiling(2.d0*pi*Rpillar/(fac*deltay)))
c     nxyp and nzp are now defined in global. Their values are obtained
c     first here

      deltaq = 2.d0*pi/dble(nxyp)
! Coordinates of tube centers
      dumx = (rlenx - 4*Rpillar)/4.d0
      dumy = (rleny - 4*Rpillar)/4.d0
      x01 = Rpillar + dumx
      y01 = Rpillar + dumy
      x02 = rlenx - x01
      y02 = rleny - y01
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l,k)
      do k=1,nzp
        do l=1,nxyp
          x_o((l-1)*nzp+k,5) = Rpillar*dcos(dble(l-1)*deltaq)+x01
          y_o((l-1)*nzp+k,5) = Rpillar*dsin(dble(l-1)*deltaq)+y01
          z_o((l-1)*nzp+k,5) = zets(k+3)

          z_o((l-1)*nzp+k,6) = zets(k+3)
          z_o((l-1)*nzp+k,8) = zets(k+3)

          x_o((l-1)*nzp+k,7) = Rpillar*dcos(dble(l-1)*deltaq)+x02
          y_o((l-1)*nzp+k,7) = Rpillar*dsin(dble(l-1)*deltaq)+y02
          z_o((l-1)*nzp+k,7) = zets(k+3)
        enddo
      enddo
C$OMP END PARALLEL DO
      n_l(1) = nx*ny
      n_l(2) = nx*ny
      n_l(3) = nx*ny
      n_l(4) = nx*ny
      n_l(5) = nxyp*nzp
      n_l(6) = nxyp*nzp
      n_l(7) = nxyp*nzp
      n_l(8) = nxyp*nzp

      num_x_o = nx ! Number of Lagrangian markers in spanwise direction
c
c      n_l_max=n_l(1)
c      do i=2,num_p
c        if (n_l(i).gt.n_l_max) n_l_max = n_l(i)
c      enddo
c      print*,'n_l_max = ',n_l_max
      print*,'n_l(1) = ',n_l(1)
      print*,'n_l(2) = ',n_l(2)
      print*,'n_l(3) = ',n_l(3)
      print*,'n_l(4) = ',n_l(4)
      print*,'n_l(5) = ',n_l(5)
      print*,'n_l(6) = ',n_l(6)
      print*,'n_l(7) = ',n_l(7)
      print*,'n_l(8) = ',n_l(8)
c      stop
c Output bed + pier
      open(1121,file="geom.dat",status='unknown')
      write(1121,*)'variables=x,y,z'
      write(1121,*)'zone t="0"',',k=',2*nx*ny+2*n_l(5)
      do i=1,nx
        do j=1,ny
          write(1121,27) x_o((i-1)*ny+j,1),y_o((i-1)*ny+j,1)
     &                  ,z_o((i-1)*ny+j,1)+rlenz/2.d0
        enddo
      enddo

      do i=1,nx
        do j=1,ny
          write(1121,27) x_o((i-1)*ny+j,3),y_o((i-1)*ny+j,3)
     &                  ,z_o((i-1)*ny+j,3)+rlenz/2.d0
        enddo
      enddo

      do l=1,nxyp
        do k=1,nzp
         write(1121,27) x_o((l-1)*nzp+k,5),y_o((l-1)*nzp+k,5)
     &                 ,z_o((l-1)*nzp+k,5)+rlenz/2.d0
        enddo
      enddo
      do l=1,nxyp
        do k=1,nzp
         write(1121,27) x_o((l-1)*nzp+k,7),y_o((l-1)*nzp+k,7)
     &                 ,z_o((l-1)*nzp+k,7)+rlenz/2.d0
        enddo
      enddo

27    format(3(3x,g25.8))
c
      return
      END SUBROUTINE initial_lag


      SUBROUTINE lagrangian_pts_gen
#include "common.inc"
c
c     This subroutine generates the Lagrangian points
c     
      INTEGER k,i,n,l,i1,i2,j1,j2,k1,k2,j
      INTEGER i_ddf,j_ddf,i_t,j_t,d
      REAL*8 ddf_dum(nx0,ny0,nz)
c      REAL*8 normal_x(nx0,ny0), normal_y(nx0,ny0), normal_z(nx0,ny0)
c      REAL*8 denom(nx,ny)

c Variables from del_fn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER zref
      REAL*8 dum1(nx),dum2(ny),dum3(nz),dum333(nz),R
      REAL*8 dzz,zminimum
      REAL*8 zarray(3),hzp,hzm,hz

c common variables
      REAL*8 FI,dum33(nz),ibmb(3),ibme(3)
c previously undeclared variables used in libraries
      REAL*8 DET
c previously undeclared variables used in libraries

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c Variables from modified_windows_fn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      INTEGER kk!,lwork,i22,i33
!      INTEGER ipiv1(n_ll)
!      PARAMETER (lwork=n_ll*(n_ll+4))
      REAL*8 ms(5),Mi(3,3)!,work(lwork)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c Variables from FI_2point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL*8 sumdum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REAL*8 eps(n_ll,num_p),ddf_control
      REAL*8 rtoto
c
c----------------------------------------------------------------------
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j) 
      do i=1,nx
        do j=1,ny
          z2d(i,j,1) = z_o((i-1)*ny+j,1) ! Bottom surface
          z2d(i,j,2) = z_o((i-1)*ny+j,3) ! Top    surface
        enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l,k)
      do k=1,nzp
        do l=1,nxyp
          rtoto = Rpillar + rf((l-1)*nzp+k,1) 
          x_o((l-1)*nzp+k,5) = rtoto*dcos(dble(l-1)*deltaq)+x01
          y_o((l-1)*nzp+k,5) = rtoto*dsin(dble(l-1)*deltaq)+y01

          if(x_o((l-1)*nzp+k,5).lt.0) then
                  x_o((l-1)*nzp+k,5) = x_o((l-1)*nzp+k,5) + rlenx
          elseif(x_o((l-1)*nzp+k,5).gt.rlenx) then
                  x_o((l-1)*nzp+k,5) = x_o((l-1)*nzp+k,5) - rlenx
          endif
          if(y_o((l-1)*nzp+k,5).lt.0) then
                  y_o((l-1)*nzp+k,5) = y_o((l-1)*nzp+k,5) + rleny
          elseif(y_o((l-1)*nzp+k,5).gt.rleny) then
                  y_o((l-1)*nzp+k,5) = y_o((l-1)*nzp+k,5) - rleny
          endif

          rtoto = Rpillar + rf((l-1)*nzp+k,2) 
          x_o((l-1)*nzp+k,7) = rtoto*dcos(dble(l-1)*deltaq)+x02
          y_o((l-1)*nzp+k,7) = rtoto*dsin(dble(l-1)*deltaq)+y02

          if(x_o((l-1)*nzp+k,7).lt.0) then
                  x_o((l-1)*nzp+k,7) = x_o((l-1)*nzp+k,7) + rlenx
          elseif(x_o((l-1)*nzp+k,7).gt.rlenx) then
                  x_o((l-1)*nzp+k,7) = x_o((l-1)*nzp+k,7) - rlenx
          endif
          if(y_o((l-1)*nzp+k,7).lt.0) then
                  y_o((l-1)*nzp+k,7) = y_o((l-1)*nzp+k,7) + rleny
          elseif(y_o((l-1)*nzp+k,7).gt.rleny) then
                  y_o((l-1)*nzp+k,7) = y_o((l-1)*nzp+k,7) - rleny
          endif
        enddo
      enddo
C$OMP END PARALLEL DO
c      
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l,k)
      do k=1,nzp
        do l=1,nxyp
          rtoto = Rpillar + rf((l-1)*nzp+k,1) + rplus
          x_o((l-1)*nzp+k,6) = rtoto*dcos(dble(l-1)*deltaq)+x01
          y_o((l-1)*nzp+k,6) = rtoto*dsin(dble(l-1)*deltaq)+y01

          if(x_o((l-1)*nzp+k,6).lt.0) then
                  x_o((l-1)*nzp+k,6) = x_o((l-1)*nzp+k,6) + rlenx
          elseif(x_o((l-1)*nzp+k,6).gt.rlenx) then
                  x_o((l-1)*nzp+k,6) = x_o((l-1)*nzp+k,6) - rlenx
          endif
          if(y_o((l-1)*nzp+k,6).lt.0) then
                  y_o((l-1)*nzp+k,6) = y_o((l-1)*nzp+k,6) + rleny
          elseif(y_o((l-1)*nzp+k,6).gt.rleny) then
                  y_o((l-1)*nzp+k,6) = y_o((l-1)*nzp+k,6) - rleny
          endif

          rtoto = Rpillar + rf((l-1)*nzp+k,2) + rplus
          x_o((l-1)*nzp+k,8) = rtoto*dcos(dble(l-1)*deltaq)+x02
          y_o((l-1)*nzp+k,8) = rtoto*dsin(dble(l-1)*deltaq)+y02

          if(x_o((l-1)*nzp+k,8).lt.0) then
                  x_o((l-1)*nzp+k,8) = x_o((l-1)*nzp+k,8) + rlenx
          elseif(x_o((l-1)*nzp+k,8).gt.rlenx) then
                  x_o((l-1)*nzp+k,8) = x_o((l-1)*nzp+k,8) - rlenx
          endif
          if(y_o((l-1)*nzp+k,8).lt.0) then
                  y_o((l-1)*nzp+k,8) = y_o((l-1)*nzp+k,8) + rleny
          elseif(y_o((l-1)*nzp+k,8).gt.rleny) then
                  y_o((l-1)*nzp+k,8) = y_o((l-1)*nzp+k,8) - rleny
          endif
        enddo
      enddo
C$OMP END PARALLEL DO
c
c
c      call ddxfd(z2d,dzdx)
c      call ddyfd(z2d,dzdy)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j) 
      do i=1,nx 
        do j=1,ny
          z_o((i-1)*ny+j,2)  = z_o((i-1)*ny+j,1) + zplus
          z_o((i-1)*ny+j,4)  = z_o((i-1)*ny+j,1) - zplus
         enddo
      enddo
C$OMP END PARALLEL DO
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j) 
c      do i=1,nx
c        do j=1,ny
c          z2d(i,j) = z_o((i-1)*ny+j,1) 
c        enddo
c      enddo
cC$OMP END PARALLEL DO
   
      do n=1,num_p
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l,i,j,k)
      do l=1,n_l(n)
        do i=1,nx0
          if (xets(i) .ge. x_o(l,n)) then
            p_iw(l,n)=i-2
            p_ie(l,n)=i+1
            goto 200
           endif
        enddo
200     do j=1,ny0
          if (yets(j) .ge. y_o(l,n)) then
            p_js(l,n)=j-2
            p_jn(l,n)=j+1
            goto 300
          endif
        enddo
300     do k=1,nz
          if (zets(k).le.z_o(l,n).and.k.ge.4.and.k.le.nz-2) then
            p_kb(l,n)=k-3
            p_kt(l,n)=k+2
            goto 400
          elseif (zets(k).le.z_o(l,n) .and. k.lt.4) then
            p_kb(l,n)=1
            p_kt(l,n)=6
            goto 400
          elseif (zets(k).le.z_o(l,n) .and. k.gt.nz-2) then
            p_kb(l,n)=nz-5
            p_kt(l,n)=nz
            goto 400
          endif
        enddo
400   enddo
C$OMP END PARALLEL DO
      enddo
c
c Initialize dums to zero
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k)
      do k=1,nz
        dum333(k)=0.d0
        dum33(k)=0.d0
        dum3(k)=0.d0
      enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i)
      do i=1,nx
        dum1(i)=0.d0
      enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j)
      do j=1,ny
        dum2(j)=0.d0
      enddo
C$OMP END PARALLEL DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c writing out the del_fn routine and its respective subroutines in
c nested do loops

c del_fn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c Setting the limits of the Lagragian points domain to reduce
c computation

      call begin_timer(100)
 
      do n=1,num_p
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,i1,i2,j,j1,j2,k,k1,k2,
C$OMP& zref,l,i_t,j_t,dum1,dum2,dum3,dum33,dum333,R,dzz,zminimum,
C$OMP& zarray,hzp,hzm,hz,sumdum,FI,Mi,ms,ibmb,DET,i_ddf,j_ddf,
C$OMP& ddf_control,ibme,kk,ddf_dum)
      do l=1,n_l(n)
        i1 = p_iw(l,n)
        i2 = p_ie(l,n)
        j1 = p_js(l,n)
        j2 = p_jn(l,n)
        k1 = p_kb(l,n)
        k2 = p_kt(l,n)
c
c Calculating the dirac delta function in the x and y directions
        do i=i1,i2
          if (i .lt. 1) then
            R=((x_o(l,n)-xets(1))
     &            +dabs(xets(nx0)-xets(nx+i)))/deltax
c            call FI_3point(R)
!%%%%%%%%%%%%%
c FI_3point(R)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (R .le. 0.5d0) then
              FI = 1.d0/3.d0*(1.d0+sqrt(-3.d0*R**2+1.d0))
            elseif (R .le. 1.5d0) then
              FI = 1.d0/6.0d0*(5.d0-3.d0*R
     &                     -sqrt(-3.d0*(1.d0-R)**2+1.d0))
            else
              FI = 0.d0
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dum1(nx+i) = FI/deltax
          elseif(i .gt. nx) then
            R=(dabs(xets(nx0)-x_o(l,n))
     &              +dabs(xets(i-nx)-xets(1)))/deltax
c            call FI_3point(R)
!%%%%%%%%%%%%%
c FI_3point(R)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (R .le. 0.5d0) then
              FI = 1.d0/3.d0*(1.d0+sqrt(-3.d0*R**2+1.d0))
            elseif (R .le. 1.5d0) then
              FI = 1.d0/6.0d0*(5.d0-3.d0*R
     &                     -sqrt(-3.d0*(1.d0-R)**2+1.d0))
            else
              FI = 0.d0
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dum1(i-nx) = FI/deltax
          else
            R=dabs((xets(i)-x_o(l,n)))/deltax
c            call FI_3point(R)
!%%%%%%%%%%%%%
c FI_3point(R)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (R .le. 0.5d0) then
              FI = 1.d0/3.d0*(1.d0+sqrt(-3.d0*R**2+1.d0))
            elseif (R .le. 1.5d0) then
              FI = 1.d0/6.0d0*(5.d0-3.d0*R
     &                     -sqrt(-3.d0*(1.d0-R)**2+1.d0))
            else
              FI = 0.d0
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dum1(i) = FI/deltax
          endif
        enddo

        do j=j1,j2
          if (j .lt. 1) then
            R=((y_o(l,n)-yets(1))
     &            +dabs(yets(ny0)-yets(ny+j)))/deltay
c            call FI_3point(R)
!%%%%%%%%%%%%%
c FI_3point(R)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (R .le. 0.5d0) then
              FI = 1.d0/3.d0*(1.d0+sqrt(-3.d0*R**2+1.d0))
            elseif (R .le. 1.5d0) then
              FI = 1.d0/6.0d0*(5.d0-3.d0*R
     &                   -sqrt(-3.d0*(1.d0-R)**2+1.d0))
            else
              FI = 0.d0
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dum2(ny+j) = FI/deltay
          elseif(j .gt. ny) then
            R=(dabs(yets(ny0)-y_o(l,n))
     &            +dabs(yets(j-ny)-yets(1)))/deltay
c            call FI_3point(R)
!%%%%%%%%%%%%%
c FI_3point(R)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (R .le. 0.5d0) then
              FI = 1.d0/3.d0*(1.d0+sqrt(-3.d0*R**2+1.d0))
            elseif (R .le. 1.5d0) then
              FI = 1.d0/6.0d0*(5.d0-3.d0*R
     &                   -sqrt(-3.d0*(1.d0-R)**2+1.d0))
            else
              FI = 0.d0
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dum2(j-ny) = FI/deltay
          else
            R=dabs((yets(j)-y_o(l,n)))/deltay
c            call FI_3point(R)
!%%%%%%%%%%%%%
c FI_3point(R)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (R .le. 0.5d0) then
              FI = 1.d0/3.d0*(1.d0+sqrt(-3.d0*R**2+1.d0))
            elseif (R .le. 1.5d0) then
              FI = 1.d0/6.0d0*(5.d0-3.d0*R
     &                   -sqrt(-3.d0*(1.d0-R)**2+1.d0))
            else
              FI = 0.d0
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dum2(j)=FI/deltay
          endif
        enddo
c
c Calculating the local dialation factor (hz) in the z direction
c ref. A. Pinelli_2010
c
c Finding the closest z cartesian point to the Lagrangian point
        zminimum = rlenz
        do k=1,nz
          dzz = dabs(zets(k)-z_o(l,n))
          if (dzz .lt. zminimum) then
            zref = k
            zminimum = dzz
          endif
        enddo
c
c defining the dialation factor hz
        if (zref .eq. 1) then
          hz=dabs(zets(zref+1)-zets(zref))/2.d0
          hzp=0.d0
          hzm=0.d0
        elseif (zref .eq. nz) then
          hz=dabs(zets(zref)-zets(zref-1))/2.d0
          hzp=0.d0
          hzm=0.d0
        elseif (zref .eq. 2) then
          zarray(1)=dabs(zets(2)-zets(1))
          zarray(2)=dabs(zets(3)-zets(2))
          zarray(3)=dabs(zets(4)-zets(3))
          hzp=max(zarray(1),zarray(2),zarray(3))
          hzm=min(zarray(1),zarray(2),zarray(3))
          hz=5.d0/6.d0*hzp+1.d0/6.d0*hzm
        else
          zarray(1) = dabs(zets(zref-1)-zets(zref-2))
          zarray(2) = dabs(zets(zref)  -zets(zref-1))
          zarray(3) = dabs(zets(zref+1)-zets(zref))
          hzp=max(zarray(1),zarray(2),zarray(3))
          hzm=min(zarray(1),zarray(2),zarray(3))
          hz=5.d0/6.d0*hzp+1.d0/6.d0*hzm
        endif
c
c Calculating the dirac delta function in the z direction
        if (zref .eq. 1 .or. zref .eq. nz) then
          do k=k1,k2
c            call FI_2point(z_o(l,n),k,zref)
!%%%%%%%%%%%%%
c FI_2point(R)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (zref .eq. 1) then
              sumdum = dabs(z_o(l,n)-zets(1))
     &                      + dabs(z_o(l,n)-zets(2))
              if (k .eq. 1) then
                FI = (zets(2)-z_o(l,n))
     &                      /sumdum/(zets(2)-zets(1))*2.d0
              elseif (k .eq. 2) then
                FI = (z_o(l,n)-zets(1))
     &                      /sumdum/(zets(3)-zets(1))*2.d0
              else
                FI = 0.d0
              endif
            endif
            if (zref .eq. nz) then
              sumdum = dabs(z_o(l,n)-zets(nz))
     &                      + dabs(z_o(l,n)-zets(nzm))
              if (k .eq. nz) then
                FI=(zets(nzm)-z_o(l,n))
     &                      /sumdum/(zets(nzm)-zets(nz))*2.d0
              elseif (k .eq. nzm) then
                FI = (z_o(l,n)-zets(nz))
     &                        /sumdum/(zets(nzmm)-zets(nz))*2.d0
              else
                FI = 0.d0
              endif
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dum3(k) = FI
          enddo
        else
          do k = k1,k2
            R=dabs(zets(k)-z_o(l,n))/hz/1.0001d0
c            call FI_3point(R)
!%%%%%%%%%%%%%
c FI_3point(R)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (R .le. 0.5d0) then
              FI = 1.d0/3.d0*(1.d0+sqrt(-3.d0*R**2+1.d0))
            elseif (R .le. 1.5d0) then
              FI = 1.d0/6.0d0*(5.d0-3.d0*R
     &                   -sqrt(-3.d0*(1.d0-R)**2+1.d0))
            else
              FI = 0.d0
            endif
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            dum3(k)=FI
          enddo
        endif
c
c Modification of dum3 for non-uniform mesh
        do k=k1,k2
          dum33(k)=dum3(k)/hz/1.0001d0
        enddo
        if (zref .eq. 1 .or. zref .eq. nz) then
          do k=k1,k2
            dum333(k)=dum3(k)
          enddo
        else
c--------------------------------------------------
cc Correction for nonuniform z direction
c            call modified_windows_fn(z_o(l,n),k1,k2,hz,dum33)
!%%%%%%%%%%%%%%%%%%%%
c modified_windows_fn
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          do k=1,5
            ms(k)=0.d0
          enddo
          do i=1,3
            ibmb(i)=0.d0
          enddo
          ibme(1)=1.d0
          ibme(2)=0.d0
          ibme(3)=0.d0
c
          do k=1,5
            do kk=k1,k2
              if (kk .eq. 1) then
                ms(k)=ms(k)+(zets(1)
     &               -z_o(l,n))**(k-1)*dum33(kk)
     &                *(zets(1)-zets(2))/2.d0
              elseif (kk .eq. nz) then
                ms(k)=ms(k)+(zets(nz)
     &               -z_o(l,n))**(k-1)*dum33(kk)
     &                *(zets(nzm)-zets(nz))/2.d0
              else
                ms(k)=ms(k)+(zets(kk)
     &               -z_o(l,n))**(k-1)*dum33(kk)
     &                *(zets(kk-1)-zets(kk+1))/2.d0
              endif
            enddo
          enddo


          do i=1,3
            do j=1,3
              Mi(i,j)=0
            enddo
          enddo
c
c          Mi(1,1)=ms(1)
c          Mi(1,2)=ms(2)
c          Mi(1,3)=ms(3)
c          Mi(2,1)=ms(2)
c          Mi(2,2)=ms(3)
c          Mi(2,3)=ms(4)
c          Mi(3,1)=ms(3)
c          Mi(3,2)=ms(4)
c          Mi(3,3)=ms(5)
cc
ccc Inverting Mi (result stored in Mi)
c          call dgetrf(3,3,Mi,3,ipiv1,i22)
c          call dgetri(3,Mi,3,ipiv1,work,lwork,i33)

ccc
! Inverting Mi (result stored in Mi)
          DET =   ms(1)*ms(3)*ms(5)
     &            - ms(1)*ms(4)*ms(4)
     &            - ms(2)*ms(2)*ms(5)
     &            + ms(2)*ms(4)*ms(3)
     &            + ms(3)*ms(2)*ms(4)
     &            - ms(3)*ms(3)*ms(3)

          Mi(1,1) = +(ms(3)*ms(5)-ms(4)*ms(4))/DET
          Mi(1,2) = -(ms(2)*ms(5)-ms(3)*ms(4))/DET
          Mi(1,3) = +(ms(2)*ms(4)-ms(3)*ms(3))/DET
          Mi(2,1) = -(ms(2)*ms(5)-ms(4)*ms(3))/DET
          Mi(2,2) = +(ms(1)*ms(5)-ms(3)*ms(3))/DET
          Mi(2,3) = -(ms(1)*ms(4)-ms(2)*ms(3))/DET
          Mi(3,1) = +(ms(2)*ms(4)-ms(3)*ms(3))/DET
          Mi(3,2) = -(ms(1)*ms(4)-ms(3)*ms(2))/DET
          Mi(3,3) = +(ms(1)*ms(3)-ms(2)*ms(2))/DET
c
cc Multiplying Mi(inverse) by ibme
          do i=1,3
            do j=1,3
              ibmb(i) = ibmb(i)+ Mi(i,j)*ibme(j)
            enddo
          enddo
cc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c--------------------------------------------------
          do k=k1,k2
            dum333(k) = (ibmb(1)+(zets(k)-z_o(l,n))*ibmb(2)
     &                     +(zets(k)-z_o(l,n))**2*ibmb(3))*dum33(k)
          enddo
        endif

        do i=i1,i2
          if(i .lt. 1) then
            i_t=nx+i
          elseif(i .gt. nx) then
            i_t=i-nx
          else
            i_t=i
          endif
          do j=j1,j2
            if(j .lt. 1) then
              j_t=ny+j
            elseif(j .gt. ny) then
              j_t=j-ny
            else
              j_t=j
            endif
            do k=k1,k2
              ddf_dum(i_t,j_t,k)=dum1(i_t)*dum2(j_t)*dum333(k)
            enddo
          enddo
        enddo
c
        do i=i1,i2
          i_ddf=i-i1+1
          if(i .lt. 1) then
            i_t=nx+i
          elseif (i .gt. nx) then
            i_t=i-nx
          else
            i_t=i
          endif
          do j=j1,j2
            j_ddf=j-j1+1
            if(j .lt.1) then
              j_t=ny+j
            elseif(j .gt. ny) then
              j_t=j-ny
            else
              j_t=j
            endif
            do k=k1,k2
              ddf(l,n,i_ddf,j_ddf,k-k1+1)=ddf_dum(i_t,j_t,k)
            enddo
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
      enddo

      call end_timer(100)
 
      do n=1,num_p,2
        do l=1,n_l(n)
          do k=1,nz
            if(zets(k) < z_o(l,n)) then
              zref = k
              goto 809
            endif
          enddo
809       continue
          eps(l,n) = deltaz(zref)
        enddo
      enddo

      do n=1,num_p,2
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l)
        do l=1,n_l(n)
         dv_l(l,n) = eps(l,n)*deltax*deltay
        enddo
C$OMP END PARALLEL DO
      enddo

c
      return
      END SUBROUTINE lagrangian_pts_gen
