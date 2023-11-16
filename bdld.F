      SUBROUTINE comp_gradient(val,dval)
#include "common.inc"

      INTEGER l
      REAL*8 val(n_ll,num_p),dval(n_ll)

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_l(1)
        dval(l)=(val(l,2) - val(l,1))/zplus
      enddo
C$OMP END PARALLEL DO

c115   format(11(1x,I8)) 
c
c      if(irkk.eq.3)then
c        write(116,1055) itime,maxval(dudn),minval(dudn)
c     *           ,sum(dudn)/dble(max(1,size(dudn))),maxval(v)
c      endif
cc      stop
c1055  format (I8,4f16.8)
c

      return
      END
c
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      SUBROUTINE comp_flux
#include "common.inc"

      INTEGER i,j,l
      REAL*8 Shields(n_ll),alpha_cr
      REAL*8 denom(nx,ny)
      REAL*8 gama(nx,ny),beta(nx,ny),alpha_x(nx,ny),alpha_y(nx,ny)
c      REAL*8 Shields_cr_1(nx,ny),Shields_cr_2(nx,ny)
c      REAL*8 dzdx_dum(nx0,ny0), dzdy_dum(nx0,ny0), dum(nx0,ny0)
      REAL*8 a1,a2
      REAL*8 slope_cr, slope_cr_b

      call begin_timer(102)

      alpha_cr = pi/6.d0
      slope_cr   = dtan(alpha_cr)
      slope_cr_b = dtan(29.d0*pi/180.d0)
       
      a1 = 4.93d0     ! Eq. 43 Wong & Parker modified MPM
      a2 = 1.6d0      ! Eq. 43 Wong & Parker modified MPM
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!           Local bed angle in 2d                 !!!!!!!!!

c Smoothen bed before computing slope 
         
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
c      do i=1,nx
c        do j=1,ny
c          dum(i,j) = z2d(i,j)
cc          dum(i,j) = dsin(3.d0*2.d0*pi*yets(j)/rleny)
cc     &             + 1.d-1*dsin(6.d1*2.d0*pi*yets(j)/rleny)
cc          z_o((i-1)*ny+j,1) = dum(i,j)
c        enddo
c      enddo
cC$OMP END PARALLEL DO
c      call gaussian(dum)

      call ddxfd(z2d,dzdx)
      call ddyfd(z2d,dzdy)

      if(mod(itime,100).eq.0.and.irkk.eq.3) then 
      print*,'z2d' ,minval(z2d(1:nx,1:ny)) + 1.d0-rlenz/2.d0
     &             ,minloc(z2d(1:nx,1:ny)) 
      print*,'z2d' ,maxval(z2d(1:nx,1:ny)) + 1.d0-rlenz/2.d0
     &             ,maxloc(z2d(1:nx,1:ny)) 
c      print*,'dum' ,minval(dum(1:nx,1:ny)) + 1.d0-rlenz/2.d0
c     &             ,minloc(dum(1:nx,1:ny)) 
c      print*,'dum' ,maxval(dum(1:nx,1:ny)) + 1.d0-rlenz/2.d0
c     &             ,maxloc(dum(1:nx,1:ny)) 
      print*,'dzdx',minval(dzdx(1:nx,1:ny))
     &             ,minloc(dzdx(1:nx,1:ny))
      print*,'dzdx',maxval(dzdx(1:nx,1:ny))
     &             ,maxloc(dzdx(1:nx,1:ny))
      print*,'dzdy',minval(dzdy(1:nx,1:ny))
     &             ,minloc(dzdy(1:nx,1:ny))
      print*,'dzdy',maxval(dzdy(1:nx,1:ny))
     &             ,maxloc(dzdy(1:nx,1:ny))
      endif

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          denom(i,j) =dsqrt(dudnx((i-1)*ny+j)**2 + dudny((i-1)*ny+j)**2)
c Compute flow orientation gama
c          dudny((i-1)*ny+j) = 5.d2
c          dudnx((i-1)*ny+j) = 0.d0
          gama(i,j) = datan2(dudnx((i-1)*ny+j),dudny((i-1)*ny+j))

c Compute bed slope x and y components 
          alpha_x(i,j) = datan(dzdx(i,j))
          alpha_y(i,j) = datan(dzdy(i,j))

c Compute bed slope along flow direction gama 
          beta(i,j) = alpha_x(i,j)*dsin(gama(i,j))  
     &              + alpha_y(i,j)*dcos(gama(i,j))

          Shields_cr_1(i,j) = Shields_cr0_1 * dcos(beta(i,j))
     &                    *(1.d0 + dtan(beta(i,j))/dtan(alpha_cr))
          Shields_cr_2(i,j) = Shields_cr0_2 * dcos(beta(i,j))
     &                    *(1.d0 + dtan(beta(i,j))/dtan(alpha_cr))
        enddo
      enddo
C$OMP END PARALLEL DO

c      open(1121,file="toto.dat",status='unknown')
c      write(1121,*)'variables=x,y,z,beta,dzdx_dum,dzdy_dum'
c      write(1121,*)'zone t="w0"',',i=',nx,',j=',ny
c      do j=1,ny
c        do i=1,nx
c         write(1121,27) deltax*dble(i-1),deltay*dble(j-1),dum(i,j) 
c     &                  + 1.d0-rlenz/2.d0,beta(i,j),dzdx_dum(i,j),
c     &                  dzdy_dum(i,j)
c        enddo
c      enddo
c      close(1121)
c27    format(6(3x,g25.8))
c
c      print*,'stop'
c       stop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                            particle 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
c      do i=1,nx
c        do j=1,ny
c          dudn((i-1)*ny+j) = dsqrt(2.d4) 
c          dudnx((i-1)*ny+j) = 1.d2
c          dudny((i-1)*ny+j) = 1.d2
c        enddo
c      enddo
cC$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          Shields((i-1)*ny+j) = c1f_1 * dudn((i-1)*ny+j)
          if(Shields((i-1)*ny+j)>Shields_cr_1(i,j)) then

            flux_x2d_1(i,j) = (dudnx((i-1)*ny+j)/denom(i,j))
     &          *(c2f_1*a1*(Shields((i-1)*ny+j)-Shields_cr_1(i,j))**a2)
     &          *phi_surf_1(i,j)/phi_cr

            flux_y2d_1(i,j) = (dudny((i-1)*ny+j)/denom(i,j))
     &          *(c2f_1*a1*(Shields((i-1)*ny+j)-Shields_cr_1(i,j))**a2)
     &          *phi_surf_1(i,j)/phi_cr

          elseif (Shields((i-1)*ny+j)<Shields_cr_1(i,j).and.
     &            dabs(dzdx(i,j))<slope_cr .and.
     &            dabs(dzdy(i,j))<slope_cr) then
            flux_x2d_1(i,j) = 0.d0
            flux_y2d_1(i,j) = 0.d0

          elseif (dabs(dzdx(i,j))>slope_cr) then 
            flux_x2d_1(i,j) =(sign(1.d0,-dzdx(i,j))/(2.d0*cfdifv(irkk)))
     &         * (0.5d0*deltax**2.d0*(dabs(dzdx(i,j))-slope_cr_b))
     &         * phi_surf_1(i,j)/phi_cr 

          elseif (dabs(dzdy(i,j))>slope_cr) then 
            flux_y2d_1(i,j) =(sign(1.d0,-dzdy(i,j))/(2.d0*cfdifv(irkk)))
     &         * (0.5d0*deltay**2.d0*(dabs(dzdy(i,j))-slope_cr_b))
     &         * phi_surf_1(i,j)/phi_cr 
          endif
        enddo
      enddo
C$OMP END PARALLEL DO
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                          particle 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          Shields((i-1)*ny+j) = c1f_2 * dudn((i-1)*ny+j) 
          if(Shields((i-1)*ny+j)>Shields_cr_2(i,j)) then
            flux_x2d_2(i,j) = (dudnx((i-1)*ny+j)/denom(i,j))
     &          *(c2f_2*a1*(Shields((i-1)*ny+j)-Shields_cr_2(i,j))**a2)
     &          *phi_surf_2(i,j)/phi_cr

            flux_y2d_2(i,j) = (dudny((i-1)*ny+j)/denom(i,j))
     &          *(c2f_2*a1*(Shields((i-1)*ny+j)-Shields_cr_2(i,j))**a2)
     &          *phi_surf_2(i,j)/phi_cr

          elseif (Shields((i-1)*ny+j)<Shields_cr_2(i,j).and.
     &            dabs(dzdx(i,j))<slope_cr .and.
     &            dabs(dzdy(i,j))<slope_cr) then
            flux_x2d_2(i,j) = 0.d0
            flux_y2d_2(i,j) = 0.d0

          elseif (dabs(dzdx(i,j))>slope_cr) then
            flux_x2d_2(i,j) =(sign(1.d0,-dzdx(i,j))/(2.d0*cfdifv(irkk)))
     &         * (0.5d0*deltax**2.d0*(dabs(dzdx(i,j))-slope_cr_b))
     &         * phi_surf_2(i,j)/phi_cr

          elseif (dabs(dzdy(i,j))>slope_cr) then
            flux_y2d_2(i,j) =(sign(1.d0,-dzdy(i,j))/(2.d0*cfdifv(irkk)))
     &         * (0.5d0*deltay**2.d0*(dabs(dzdy(i,j))-slope_cr_b))
     &         * phi_surf_2(i,j)/phi_cr

          endif
        enddo
      enddo
C$OMP END PARALLEL DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                             Total Flux 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          flux_x2d(i,j) = flux_x2d_1(i,j) + flux_x2d_2(i,j)  
          flux_y2d(i,j) = flux_y2d_1(i,j) + flux_y2d_2(i,j)  
          flux_2d(i,j)  = dsqrt(flux_x2d(i,j)**2 + flux_y2d(i,j)**2)
        enddo
      enddo
C$OMP END PARALLEL DO


c
      call end_timer(102)
c

      return
      END
c
      SUBROUTINE extract_surf_volfrac 
#include "common.inc"

      INTEGER i,j,k,toto 

      if (itime.eq.istart+1)then
        do i=1,nx
          do j=1,ny
            do k=1,nz_phi
              if(z2d(i,j).lt.zets_phi(k))then 
                k_bed(i,j) = k-1
                toto       = k_bed(i,j) 
                phi_surf_1(i,j) = phi_1(i,j,k_bed(i,j))
                phi_surf_2(i,j) = phi_2(i,j,k_bed(i,j))
                if(phi_surf_1(i,j)+phi_surf_2(i,j).eq.0) then
711               toto = toto - 1
                  phi_surf_1(i,j) = phi_1(i,j,toto)
                  phi_surf_2(i,j) = phi_2(i,j,toto)
c                  if(phi_surf_1(i,j)+phi_surf_2(i,j).eq.0) goto 711 
                  if(phi_surf_1(i,j)+phi_surf_2(i,j).eq.0) print*,'Stop'
                endif
                goto 809
              endif
            enddo
809       continue
          enddo
        enddo

      else !use value of bed from previous time step 

        do i=1,nx
          do j=1,ny
            do k=k_bed(i,j),k_bed(i,j)+2
              if(z2d(i,j).lt.zets_phi(k))then
                k_bed(i,j) = k-1
                toto       = k_bed(i,j) 
                phi_surf_1(i,j) = phi_1(i,j,k_bed(i,j))
                phi_surf_2(i,j) = phi_2(i,j,k_bed(i,j))
                if(phi_surf_1(i,j)+phi_surf_2(i,j).eq.0) then
712               toto = toto - 1
                  phi_surf_1(i,j) = phi_1(i,j,toto)
                  phi_surf_2(i,j) = phi_2(i,j,toto)
c                  if(phi_surf_1(i,j)+phi_surf_2(i,j).eq.0) goto 712
                  if(phi_surf_1(i,j)+phi_surf_2(i,j).eq.0) print*,'Stop'
                endif
                goto 810
              endif
            enddo
810       continue
          enddo
        enddo
      endif

      return
      END

      SUBROUTINE comp_height
#include "common.inc"

      INTEGER i,j,l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!       Update one spanwise x-plane at a time     !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call begin_timer(103)

      call advance_h 
c      call sandslide
      call conserve_h 
      call extract_surf_volfrac
      call update_phi 

      call end_timer(103)
c
      return
      END
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine conserve_h
#include "common.inc"
      INTEGER l,i,j
      Real*8 h_sum,h2d(nx,ny)

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          if (z_o((i-1)*ny+j,1).lt.(-rlenz/2.d0 + 5.d-3)) then
            z_o((i-1)*ny+j,1) = -rlenz/2.d0+5.d-3
          endif
        enddo
      enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          h2d(i,j)=z_o((i-1)*ny+j,1)
        enddo
      enddo
C$OMP END PARALLEL DO
      h_sum = sum(h2d)
c      print*,'h_sum = ',h_sum

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          z_o((i-1)*ny+j,1) = z_o((i-1)*ny+j,1)
     &           -(h_sum-h0_sum)/dble(nx*ny) 
        enddo
      enddo
C$OMP END PARALLEL DO


      return
      END
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine update_phi 
#include "common.inc"
      real *8 dqxdx_i( nx0,ny0), dqydy_i( nx0,ny0)
      real *8 exrc(3), exzc(3) 
      real *8 toto_1(nx0,ny0), toto_2(nx0,ny0)
      real *8 sum_1, sum_2
      INTEGER i,j,n

      exrc(1) =  0.d0
      exrc(2) = -5.d0/9.d0
      exrc(3) = -153.d0/128.d0

      exzc(1) = 1.d0/3.d0
      exzc(2) = 15.d0/16.d0
      exzc(3) = 8.d0/15.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!     Compute bed & flux derivatives     !!!!!!!!!!!!!!
!!!!!!!!!!!!!!             particle 1                 !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call ddxfd(flux_x2d_1,dqxdx_i)
      call ddyfd(flux_y2d_1,dqydy_i)

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do i=1,nx
          do j=1,ny
           exner_rhs_1(i,j) = exrc(irkk)*exner_rhs_1(i,j)
     &        - (dt/deltaz_phi)
     &        * (dqxdx_i(i,j)-(phi_surf_1(i,j)/phi_cr)*dqxdx(i,j)
     &         + dqydy_i(i,j)-(phi_surf_1(i,j)/phi_cr)*dqydy(i,j))

           phi_surf_1(i,j)=phi_surf_1(i,j)+exzc(irkk)*exner_rhs_1(i,j)
           if(phi_surf_1(i,j).lt.0.d0)   phi_surf_1(i,j) = 0.d0
           if(phi_surf_1(i,j).gt.phi_cr) phi_surf_1(i,j) = phi_cr
          toto_1(i,j) = phi_surf_1(i,j)
          enddo
        enddo
C$OMP END PARALLEL DO
c
        ! nadim
      if(mod(itime,100).eq.0.and.irkk.eq.3) then 
        print*,itime,irkk,maxval(exner_rhs_1(1:nx,1:ny))
     &                   ,maxloc(exner_rhs_1(1:nx,1:ny))
        print*,itime,irkk,minval(exner_rhs_1(1:nx,1:ny))
     &                   ,minloc(exner_rhs_1(1:nx,1:ny))
      endif
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!     Compute bed & flux derivatives     !!!!!!!!!!!!!!
!!!!!!!!!!!!!!             particle 2                 !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call ddxfd(flux_x2d_2,dqxdx_i)
      call ddyfd(flux_y2d_2,dqydy_i)

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do i=1,nx
          do j=1,ny
           exner_rhs_2(i,j) = exrc(irkk)*exner_rhs_2(i,j)
     &        - (dt/deltaz_phi)
     &        * (dqxdx_i(i,j)-(phi_surf_2(i,j)/phi_cr)*dqxdx(i,j)
     &         + dqydy_i(i,j)-(phi_surf_2(i,j)/phi_cr)*dqydy(i,j))

           phi_surf_2(i,j)=phi_surf_2(i,j)+exzc(irkk)*exner_rhs_2(i,j)
           if(phi_surf_2(i,j).lt.0.d0)   phi_surf_2(i,j) = 0.d0
           if(phi_surf_2(i,j).gt.phi_cr) phi_surf_2(i,j) = phi_cr
          toto_2(i,j) = phi_surf_2(i,j)
          enddo
        enddo
C$OMP END PARALLEL DO
      if(mod(itime,100).eq.0.and.irkk.eq.3) then 
        print*,itime,irkk,maxval(exner_rhs_2(1:nx,1:ny))
     &                   ,maxloc(exner_rhs_2(1:nx,1:ny))
        print*,itime,irkk,minval(exner_rhs_2(1:nx,1:ny))
     &                   ,minloc(exner_rhs_2(1:nx,1:ny))

        print*, maxval(phi_surf_1(1:nx,1:ny))
        print*, minval(phi_surf_2(1:nx,1:ny))
        print*, maxval(phi_surf_2(1:nx,1:ny))
        print*, minval(phi_surf_1(1:nx,1:ny))
      endif
c        if(irkk.eq.3) stop
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do i=1,nx
          do j=1,ny
            phi_surf_1(i,j) =phi_surf_1(i,j)
     &                       *phi_cr/(toto_1(i,j)+toto_2(i,j)) 
c            phi_1(i,j,k_bed(i,j)) = phi_surf_1(i,j)
c            phi_1(i,j,k_bed(i,j)-1) = 0.d0
c            phi_1(i,j,k_bed(i,j)-2) = 0.d0

            phi_surf_2(i,j) =phi_surf_2(i,j)
     &                       *phi_cr/(toto_1(i,j)+toto_2(i,j)) 
c            phi_2(i,j,k_bed(i,j)) = phi_surf_2(i,j)
c            phi_2(i,j,k_bed(i,j)-1) = 0.d0
c            phi_2(i,j,k_bed(i,j)-2) = 0.d0
          enddo
        enddo
C$OMP END PARALLEL DO
c
c        print*, maxval(phi_surf_1(1:nx,1:ny))
c        print*, minval(phi_surf_2(1:nx,1:ny))
c        print*, maxval(phi_surf_2(1:nx,1:ny))
c        print*, minval(phi_surf_1(1:nx,1:ny))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c                       Conserve phi_surf_1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
        do n=1,3
c        if( mod((itime-istart),oi_mean/2).eq.0) then
          sum_1 = sum(phi_surf_1(1:nx,1:ny))
          sum_2 = sum(phi_surf_2(1:nx,1:ny))
c          sum_t = phi_cr*dble(nx*ny)

c        print*,k_bed0, k_bed(1,1)
c        print*,phi_cr
c        stop

c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
          do i=1,nx
            do j=1,ny
              phi_surf_1(i,j)=phi_surf_1(i,j)*(frac_1*sum_t/sum_1)
              toto_1(i,j) = phi_surf_1(i,j)
  
              phi_surf_2(i,j)=phi_surf_2(i,j)*(frac_2*sum_t/sum_1)
              toto_2(i,j) = phi_surf_2(i,j)
            enddo
          enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
          do i=1,nx
            do j=1,ny
              phi_surf_1(i,j) =phi_surf_1(i,j)
     &                       *phi_cr/(toto_1(i,j)+toto_2(i,j)) 
              phi_surf_2(i,j) =phi_surf_2(i,j)
     &                       *phi_cr/(toto_1(i,j)+toto_2(i,j)) 
            enddo
          enddo
C$OMP END PARALLEL DO
cc
        enddo ! n=1,3
c
c        print*, maxval(phi_surf_1(1:nx,1:ny))
c        print*, minval(phi_surf_2(1:nx,1:ny))
c        print*, maxval(phi_surf_2(1:nx,1:ny))
c        print*, minval(phi_surf_1(1:nx,1:ny))
cc
        sum_1 = sum(phi_surf_1(1:nx,1:ny))
        sum_2 = sum(phi_surf_2(1:nx,1:ny))
c
      if(mod(itime,100).eq.0.and.irkk.eq.3) then 
        print*,'itime = ',itime,'sum_1 = ',sum_1
        print*,'itime = ',itime,'sum_2 = ',sum_2
      endif
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
          do i=1,nx
            do j=1,ny
              phi_1(i,j,k_bed(i,j)) = phi_surf_1(i,j)
              phi_1(i,j,k_bed(i,j)+1) = 0.d0
              phi_1(i,j,k_bed(i,j)+2) = 0.d0
              phi_1(i,j,k_bed(i,j)+3) = 0.d0
  
              phi_2(i,j,k_bed(i,j)) = phi_surf_2(i,j)
              phi_2(i,j,k_bed(i,j)+1) = 0.d0
              phi_2(i,j,k_bed(i,j)+2) = 0.d0
              phi_2(i,j,k_bed(i,j)+3) = 0.d0
            enddo
          enddo
C$OMP END PARALLEL DO
c        endif ! oi_mean
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      return
      END
c
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine advance_h
#include "common.inc"
      real *8 exrc(3), exzc(3) 
      INTEGER i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!  Compute the mean flux    !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      q_mean = sum(dabs(flux_2d))/dble(nx*ny)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!     Compute bed & flux derivatives     !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call ddxfd(flux_x2d,dqxdx)
      call ddyfd(flux_y2d,dqydy)
c      call ddxfd(z2d,dzdx)
c      call ddyfd(z2d,dzdy)
      call d2dx2fd(z2d,d2zdx2)
      call d2dy2fd(z2d,d2zdy2)
c      call ddx2d(dzdx,d2zdx2)
c      call ddy2d(dzdy,d2zdy2)

      exrc(1) =  0.d0
      exrc(2) = -5.d0/9.d0
      exrc(3) = -153.d0/128.d0

      exzc(1) = 1.d0/3.d0
      exzc(2) = 15.d0/16.d0
      exzc(3) = 8.d0/15.d0

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do i=1,nx
          do j=1,ny
            exner_rhs(i,j) = exrc(irkk)*exner_rhs(i,j)-(dt/phi_cr)
     &        * ((dqxdx(i,j)+dqydy(i,j))
     &        - eps_bed*q_mean*(d2zdx2(i,j)+d2zdy2(i,j)))
            z2d(i,j) = z2d(i,j) + exzc(irkk)*exner_rhs(i,j)
            z_o((i-1)*ny+j,1) = z2d(i,j)
          enddo
        enddo
C$OMP END PARALLEL DO
c
      return
      END
c
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine sandslide
#include "common.inc"

      integer n,i,j,ip1,im1,jp1,jm1
      REAL*8 slope_cr1, slope_cr2

      call ddxfd(z2d,dzdx)
      call ddyfd(z2d,dzdy)

*======================================================================
      slope_cr1 = dtan(30.d0*pi/180.d0)
      slope_cr2 = dtan(29.d0*pi/180.d0)

      n=0
      do while(max(maxval(abs(dzdx)),maxval(abs(dzdy))).ge.slope_cr1
     &             .and.n.le.5)
*======================================================================
***************** Sand-slide downstream direction *********************
*======================================================================
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jp1)
        do i=1,nx
          do j=1,ny
            jp1=j+1
            if(j.eq.ny)then
              jp1=1
            endif
            if(z2d(i,j) > z2d(i,jp1) + slope_cr1*deltay) then
              z2d(i,j) = 0.5d0*(z2d(i,j)+z2d(i,jp1) + slope_cr2*deltay)
              z2d(i,jp1) = z2d(i,j) - slope_cr2*deltay
            endif
          enddo
        enddo
C$OMP END PARALLEL DO
c
*======================================================================
***************** Sand-slide cross-stream + direction *****************
*======================================================================
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,ip1)
        do j=1,ny
          do i=1,nx
            ip1=i+1
            if(i.eq.nx)then
              ip1=1
            endif
            if(z2d(i,j) > z2d(ip1,j) + slope_cr1*deltax) then
              z2d(i,j) = 0.5d0*(z2d(i,j)+z2d(ip1,j) + slope_cr2*deltax)
              z2d(ip1,j) = z2d(i,j) - slope_cr2*deltax
            endif
          enddo
        enddo
C$OMP END PARALLEL DO
c
*======================================================================
****************** Sand-slide upstream direction ********************
*======================================================================
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,jm1)
        do i=1,nx
          do j=1,ny
            jm1=j-1
            if(j.eq.1)then
              jm1=ny
            endif
            if(z2d(i,j) > z2d(i,jm1) + slope_cr1*deltay) then
              z2d(i,jm1) = 0.5d0*(z2d(i,j)+z2d(i,jm1)-slope_cr2*deltay)
              z2d(i,j) = z2d(i,jm1) + slope_cr2*deltay
            endif
          enddo
        enddo
C$OMP END PARALLEL DO
*======================================================================
***************** Sand-slide cross-stream - direction *****************
*======================================================================
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,im1)
        do j=1,ny
          do i=1,nx
            im1=i-1
            if(i.eq.1)then
              im1=nx
            endif
            if(z2d(i,j) > z2d(im1,j) + slope_cr1*deltax) then
              z2d(im1,j) = 0.5d0*(z2d(i,j)+z2d(im1,j)-slope_cr2*deltax)
              z2d(i,j) = z2d(im1,j) + slope_cr2*deltax
            endif
          enddo
        enddo
C$OMP END PARALLEL DO
c
      n=n+1
      call ddxfd(z2d,dzdx)
      call ddyfd(z2d,dzdy)
c      print*,itime,n,max(maxval(abs(dzdx)),maxval(abs(dzdy)))
      enddo !do while
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
          do i=1,nx
            z_o((i-1)*ny+j,1) = z2d(i,j)
          enddo
        enddo
C$OMP END PARALLEL DO
c
      return
      END
c
      subroutine gaussian(dum_new)
#include "common.inc"
      INTEGER n,i,im2,im1,ip1,ip2,j,jm2,jm1,jp1,jp2
      Real*8 dum(nx0,ny0),dum_new(nx0,ny0)
      INTEGER n_max

c      print*,'maxval(dum) = ',maxval(dum_new(1:nx,1:ny))
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do i=1,nx
        do j=1,ny
          dum(i,j) =z_o((i-1)*ny+j,1) 
        enddo
      enddo
C$OMP END PARALLEL DO

      n_max = 1
      do n=1,n_max
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,im2,im1,ip1,ip2,j,jm2,jm1,
C$OMP& jp1,jp2)
        do i=1,nx
          im2=i-2
          im1=i-1
          ip1=i+1
          ip2=i+2
          if(i.eq.1)then
            im1=nx
            im2=nx-1
          elseif(i.eq.2)then
            im2=nx
          elseif(i.eq.nx)then
            ip1=1
            ip2=2
          elseif(i.eq.nx-1)then
            ip2=1
          endif
          do j=1,ny
            jm2=j-2
            jm1=j-1
            jp1=j+1
            jp2=j+2
            if(j.eq.1)then
              jm1=ny
              jm2=ny-1
            elseif(j.eq.2)then
              jm2=ny
            elseif(j.eq.ny)then
              jp1=1
              jp2=2
            elseif(j.eq.ny-1)then
              jp2=1
            endif
            dum_new(i,j) = (1./273.)*(41.*dum(i,j) + 26.*
     &     (dum(ip1,j)  +dum(im1,j)  +dum(i,jp1)  +dum(i,jm1))  +16.*
     &     (dum(ip1,jp1)+dum(ip1,jm1)+dum(im1,jp1)+dum(im1,jm1))+ 7.*
     &     (dum(ip2,j)  +dum(im2,j)  +dum(i,jp2)  +dum(i,jm2))  + 4.*
     &     (dum(ip2,jp1)+dum(ip2,jm1)+dum(ip1,jp2)+dum(ip1,jm2))+ 4.*
     &     (dum(im2,jp1)+dum(im2,jm1)+dum(im1,jp2)+dum(im1,jm2))+ 1.*
     &     (dum(ip2,jp2)+dum(ip2,jm2)+dum(im2,jp2)+dum(im2,jm2)) )
          enddo ! ny
        enddo ! nx
C$OMP END PARALLEL DO

        if (n<n_max) then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
          do i=1,nx
            do j=1,ny
            dum(i,j) = dum_new(i,j)
            enddo
          enddo
C$OMP END PARALLEL DO
        endif !(n<n_max)
      enddo ! do n
c
      return
      END
c
*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c

