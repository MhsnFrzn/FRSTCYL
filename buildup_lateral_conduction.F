      SUBROUTINE comp_gradient(val,dval)
#include "common.inc"

      INTEGER l
      REAL*8 val(n_ll,num_p),dval(n_ll)

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        dval(l)=(val(l,2) - val(l,1))/zplus
      enddo
C$OMP END PARALLEL DO

      return
      END SUBROUTINE comp_gradient
cc
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
cc
      SUBROUTINE fro_boundary(val1,val2)
#include "common.inc"

      INTEGER l
      REAL*8 val1(n_ll), val2(n_ll), T_fro(n_ll), rho_air(n_ll)
      REAL*8 k_fro(n_ll), Latent(n_ll), massdiff(n_ll), p_fro(n_ll)
      REAL*8 k_air(n_ll), eps(n_ll), T_cel(n_ll), rho_fro(n_ll)
      REAL*8 rho_fro_n(n_ll), T_fro_n(n_ll)
      REAL*8 const1(n_ll), const2(n_ll), const3(n_ll), const4(n_ll)
      REAL*8 var1(n_ll), var2(n_ll), var3(n_ll)

      REAL*8 ck1, ck2, ck3, ck4 
      REAL*8 ce0, ce1, ce2, ce3, ce4, ce5
      REAL*8 cp0, cp1, cp2, cp3, cp4, cp5, cp6
      REAL*8 p_atm, R_sgc, nu, cp

c Coefficients to determine air density 
      p_atm = 101.325d3 ! N/m^2
      R_sgc = 287.058d0 ! specific gas constant J/(kg.K) 
      nu        = 1.42d-5! kinematic viscosity of air

c Coefficients to determine air conductivity 
      ck1 =  1.5207d-11
      ck2 = -4.8574d-8
      ck3 =  1.0184d-4
      ck4 = -3.9333d-4

c Coefficients to determine eps from Mago & Sherif 
      ce0 =  1.00391d0 
      ce1 = -7.82205d-6
      ce2 =  6.94682d-7
      ce3 =  3.04059d-9
      ce4 = -2.78520d-11
      ce5 = -5.65600d-13

c Coefficients to determine p_fro from Mago & Sherif 
      cp0 = -5.6745359d3 
      cp1 =  6.3925247d0 
      cp2 = -0.9677843d-2 
      cp3 =  6.2215701d-7 
      cp4 =  2.0747825d-9
      cp5 = -9.4840240d-13
      cp6 =  4.1635019d0


      cp        = 2.04d3 ! specific heat of ice
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c Variables evaluated at old time step
        T_fro(l)   = T_wall + (T_inf-T_wall)*tt_p(l,1)!Frost surface temperature
c        T_fro(l) = max(T_wall,T_fro(l)) 
c Temperature cannot go below wall temperature (Case when Nu<0)
        rho_air(l) = p_atm/(R_sgc*T_fro(l)) ! rho=p/(R*T)
        k_air(l)   = ck1*T_fro(l)**3+ck2*T_fro(l)**2+ck3*T_fro(l)+ck4
        rho_fro(l) = 650.d0*dexp(0.227d0*(T_fro(l)-273.15d0)) !Mago & Sherif 
        k_fro(l)   = 2.4248d-2 + 7.2311d-4*rho_fro(l) 
     &                         + 1.1830d-6*rho_fro(l)**2
        Latent(l)  = -0.1083d0*(1.8d0*(T_fro(l)-273.15d0)+32.d0)+2833.d0
        massdiff(l)= 0.23d-4 ! We could find a fit for this as well



c        tt_p(l,1) = ((z_o(l,1)+z_o_avg-rlenz/2.d0+1.d0) 
c     &               * (k_air(l)/k_fro(l)) * val1(l) + 
c     &               ((z_o(l,1)+z_o_avg-rlenz/2.d0+1.d0) * Y_inf 
c     &               * rho_air(l) * Latent(l) * massdiff(l) / (k_fro(l)
c     &               * (T_inf-T_wall))) * val2(l)) 

      enddo
C$OMP END PARALLEL DO
c        tt_p(l,1) = max(tt_p(l,1),0.d0)
      if (irkk.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,n_ll
          var1(l) = cp*(z_o(l,1)+z_o_avg-rlenz/2.d0+1.d0)*rho_fro(l)
     &                * re*nu/2.d0
          var2(l) = k_air(l)*val1(l)  
     &                + val2(l)*Y_inf*rho_air(l)*Latent(l)*massdiff(l)/
     &                  (T_inf-T_wall)
          var3(l) = k_fro(l)/(z_o(l,1)+z_o_avg-rlenz/2.d0+1.d0)

          tt_p(l,1) = (var1(l)*tt_p(l,1) + var2(l)*dt)
     &               /(var1(l)           + var3(l)*dt)
        enddo
C$OMP END PARALLEL DO
      endif
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c Compute temperature and density at new time step
        T_fro_n(l)   = T_wall + (T_inf-T_wall)*tt_p(l,1)!Air temp at frost layer 
c        T_fro_n(l) = max(T_wall,T_fro_n(l)) 
c Temperature cannot go below wall temperature (Case when Nu<0)
        rho_fro_n(l) = 650.d0*dexp(0.227d0*(T_fro_n(l)-273.15d0)) !Mago & Sherif 
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
        p_fro(l) = dexp(cp0/T_fro(l) +cp1 +cp2*T_fro(l) +cp3*T_fro(l)**2
     &                 +cp4*T_fro(l)**3 +cp5*T_fro(l)**4 
     &                 +cp6*dlog(T_fro(l))) 

        T_cel(l) = T_fro(l) - 273.15d0
        eps(l) = ce0 +ce1*T_cel(l)    + ce2*T_cel(l)**2 +ce3*T_cel(l)**3
     &               +ce4*T_cel(l)**4 + ce5*T_cel(l)**5 

        tt2_p(l,1) = (eps(l) * (p_fro(l)/p_atm) / Y_inf)
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      enddo
C$OMP END PARALLEL DO

      if(mod(itime,oi_mean).eq.0.and.irkk.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
         do l=1,n_ll
           const1(l) = (z_o(l,1)+z_o_avg+1.d0-rlenz/2.d0)
     &       *k_air(l)/k_fro(l)
           const2(l) = (z_o(l,1)+z_o_avg+1.d0-rlenz/2.d0)*Y_inf
     &       *rho_air(l)*Latent(l)*massdiff(l)/(k_fro(l)*(T_inf-T_wall))
           const3(l) = (rho_air(l)/rho_fro(l))*Y_inf/(re*pr2)
           const4(l) = rho_fro(l)/rho_fro_n(l)
         enddo 
C$OMP END PARALLEL DO
         print*,'max const1',maxval(const1)
         print*,'min const1',minval(const1)
         print*,'max const2',maxval(const2)
         print*,'min const2',minval(const2)
         print*,'max const3',maxval(const3)
         print*,'min const3',minval(const3)
         print*,'max const4',maxval(const4)
         print*,'min const4',minval(const4)
 
         print*,'max tt_p',maxval(tt_p(:,1))
         print*,'min tt_p',minval(tt_p(:,1))
         print*,'avg tt_p',sum(tt_p(:,1))/dble(n_ll)
         print*,'max tt2_p',maxval(tt2_p(:,1))
         print*,'min tt2_p',minval(tt2_p(:,1))
       endif
c        print*,'max T_fro',maxval(T_fro)
c        print*,'max rho_air',maxval(rho_air)
c        print*,'max k_air',maxval(k_air)
c        print*,'max rho_fro',maxval(rho_fro)
c        print*,'max rho_fro_n',maxval(rho_fro_n)
c        print*,'max k_fro',maxval(k_fro)
c        print*,'max Latent',maxval(Latent)
c        print*,'max massdiff',maxval(massdiff)
c
c        print*,'min T_fro',minval(T_fro)
c        print*,'min rho_air',minval(rho_air)
c        print*,'min k_air',minval(k_air)
c        print*,'min rho_fro',minval(rho_fro)
c        print*,'min rho_fro_n',minval(rho_fro_n)
c        print*,'min k_fro',minval(k_fro)
c        print*,'min Latent',minval(Latent)
c        print*,'min massdiff',minval(massdiff)


c        if(itime.eq.3) stop

      call update_h(rho_air,rho_fro,rho_fro_n,val2)
      call lagrangian_pts_gen
c
      return
      END SUBROUTINE fro_boundary
cc
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine update_temp(val1,rho_fro,T_fro,k_air,k_fro)
#include "common.inc"
      REAL*8 exrc(3), exzc(3) 
      REAL*8 k_air(n_ll), rho_fro(n_ll), val1(n_ll),k_fro(n_ll)
      REAL*8 rho_layer(n_ll),T_fro(n_ll)
      REAL*8 z_o_c(n_ll)
      REAL*8 temp_2d(nx0,ny0), d2tdx2(nx0,ny0), d2tdy2(nx0,ny0)
      REAL*8 diff_coef,nu,cp,eps
      INTEGER l,i,j

      diff_coef = 1.d-2
      nu        = 1.42d-5! kinematic viscosity of air
      cp        = 2.04d3 ! specific heat of ice
      eps       = 1.d-6  ! small number 

      exrc(1) =  0.d0
      exrc(2) = -5.d0/9.d0
      exrc(3) = -153.d0/128.d0

      exzc(1) = 1.d0/3.d0
      exzc(2) = 15.d0/16.d0
      exzc(3) = 8.d0/15.d0


C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        rho_layer(l) = 0.5d0*(650.d0*dexp(0.227d0*(T_wall-273.15d0))
     &                       +650.d0*dexp(0.227d0*(T_fro(l)-273.15d0)))
        z_o_c(l)     = z_o(l,1) + z_o_avg + 1.d0-rlenz/2.d0 + eps
      enddo
C$OMP END PARALLEL DO

c 1 index to 2 index representation
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j=1,ny
        do i=1,nx
          temp_2d(i,j) = tt_p((i-1)*ny+j,1)
        enddo
      enddo
C$OMP END PARALLEL DO

      call d2dx2fd(temp_2d,d2tdx2)
      call d2dy2fd(temp_2d,d2tdy2)

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        temp_rhs(l) = exrc(irkk)*temp_rhs(l) 
     &     +AF*dt*((4.d0/(cp*z_o_c(l)*rho_layer(l)*re*nu))
     &      *(k_air(l)*(1.d0-tt_p(l,1))*val1(l)
     &          -k_fro(l)*tt_p(l,1)/(z_o_c(l))) 
     &     +diff_coef*d2tdx2(1+int((l-1)/ny),1+mod(l-1,ny))
     &     +diff_coef*d2tdy2(1+int((l-1)/ny),1+mod(l-1,ny)))

        tt_p(l,1) =  tt_p(l,1) + exzc(irkk)*temp_rhs(l)

c        z_o_c(l) = z_o_c(l)-(1.d0-rho_fro(l)/rho_fro_n(l)) 
c     &           * (z2d(1+int((l-1)/ny),1+mod(l-1,ny))+ z_o_avg
c     &           + 1.d0-rlenz/2.d0)
      enddo
C$OMP END PARALLEL DO

      print*,itime,irkk
      print*,maxval(tt_p(:,1))
      print*,minval(tt_p(:,1))
c

c        print*,itime,irkk
c        print*,'z_o_avg',z_o_avg
c        print*,'min z_o_c',minval(z_o_c)  + 1.d0-rlenz/2.d0
c        print*,'max z_o_c',maxval(z_o_c)  + 1.d0-rlenz/2.d0
c        print*,'min z_o',minval(z_o(:,1)) + 1.d0-rlenz/2.d0
c        print*,'max z_o',maxval(z_o(:,1)) + 1.d0-rlenz/2.d0
c        print*,'avg z_o',sum(z_o(:,1))/dble(n_ll) + 1.d0-rlenz/2.d0
c        print*,'avg z_o_c',sum(z_o_c)/dble(n_ll) + 1.d0-rlenz/2.d0
c        print*,'delta z_o',maxval(z_o_c)-minval(z_o_c)
c        if(itime.eq.3) stop
c
      return
      END SUBROUTINE update_temp
cc
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine update_h(rho_air,rho_fro,rho_fro_n,val2)
#include "common.inc"
      REAL*8 exrc(3), exzc(3) 
      REAL*8 rho_air(n_ll), rho_fro(n_ll), val2(n_ll),rho_fro_n(n_ll)
      REAL*8 z_o_c(n_ll),z_o_c_old(n_ll)
c      REAL*8 z_o_2d(nx0,ny0), d2zdx2(nx0,ny0), d2zdy2(nx0,ny0)
c      REAL*8 diff_coef
      INTEGER l!,i,j

c      diff_coef = 1.d-2

      exrc(1) =  0.d0
      exrc(2) = -5.d0/9.d0
      exrc(3) = -153.d0/128.d0

      exzc(1) = 1.d0/3.d0
      exzc(2) = 15.d0/16.d0
      exzc(3) = 8.d0/15.d0


C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        z_o_c(l)     = z_o(l,1) + z_o_avg
        z_o_c_old(l) = z_o_c(l) 
      enddo
C$OMP END PARALLEL DO

cc 1 index to 2 index representation
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
c      do j=1,ny
c        do i=1,nx
c          z_o_2d(i,j) = z_o((i-1)*ny+j,1)
c        enddo
c      enddo
cC$OMP END PARALLEL DO
c    
c      call d2dx2fd(z_o_2d,d2zdx2)
c      call d2dy2fd(z_o_2d,d2zdy2)
     

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        exner_rhs(l) = exrc(irkk)*exner_rhs(l) 
     &     +AF*dt*(rho_air(l)/rho_fro(l)) * (Y_inf/(re*pr2)) * val2(l)
c     &     +diff_coef*d2zdx2(1+int((l-1)/ny),1+mod(l-1,ny))
c     &     +diff_coef*d2zdy2(1+int((l-1)/ny),1+mod(l-1,ny)))

        z_o_c(l) =  z_o_c(l) + exzc(irkk)*exner_rhs(l)

c        z_o_c(l) = z_o_c(l)-(1.d0
c     &           -min(rho_fro(l)/rho_fro_n(l),1.d0))
c     &           *(z_o_c_old(l) + 1.d0-rlenz/2.d0)

c        if(z_o_c(l).lt.z_o_c_old(l)) z_o_c(l) = z_o_c_old(l)

c        z_o_c(l) = z_o_c(l)-(1.d0-rho_fro(l)/rho_fro_n(l)) 
c     &           * (z2d(1+int((l-1)/ny),1+mod(l-1,ny))+ z_o_avg
c     &           + 1.d0-rlenz/2.d0)
      enddo
C$OMP END PARALLEL DO

c Compaction
      if(irkk.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,n_ll
          z_o_c(l) = z_o_c(l)-(1.d0-rho_fro(l)/rho_fro_n(l))
     &           *(z_o_c_old(l) + 1.d0-rlenz/2.d0)
        enddo
C$OMP END PARALLEL DO
      endif
c
c      if(irkk.eq.3) print*,z_o_c(1),z_o_c_old(1)
c      if(irkk.eq.3) print*,z_o_c(1) + 1.d0-rlenz/2.d0,
c     &                 z_o_c_old(1) + 1.d0-rlenz/2.d0
c      if(irkk.eq.3) stop

      z_o_avg = sum(z_o_c)/dble(n_ll)+1.d0-rlenz/2.d0


C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        z_o(l,1) = z_o_c(l) - z_o_avg
      enddo
C$OMP END PARALLEL DO


c        print*,itime,irkk
c        print*,'z_o_avg',z_o_avg
c        print*,'min z_o_c',minval(z_o_c)  + 1.d0-rlenz/2.d0
c        print*,'max z_o_c',maxval(z_o_c)  + 1.d0-rlenz/2.d0
c        print*,'min z_o',minval(z_o(:,1)) + 1.d0-rlenz/2.d0
c        print*,'max z_o',maxval(z_o(:,1)) + 1.d0-rlenz/2.d0
c        print*,'avg z_o',sum(z_o(:,1))/dble(n_ll) + 1.d0-rlenz/2.d0
c        print*,'avg z_o_c',sum(z_o_c)/dble(n_ll) + 1.d0-rlenz/2.d0
c        print*,'delta z_o',maxval(z_o_c)-minval(z_o_c)
c        if(itime.eq.3) stop



      if(mod(itime,oi_mean).eq.0.and.irkk.eq.3) then
        print*,'z_o_avg',z_o_avg
        print*,'min z_o_c',minval(z_o_c)  + 1.d0-rlenz/2.d0
        print*,'max z_o_c',maxval(z_o_c)  + 1.d0-rlenz/2.d0
        print*,'min z_o',minval(z_o(:,1)) + 1.d0-rlenz/2.d0
        print*,'max z_o',maxval(z_o(:,1)) + 1.d0-rlenz/2.d0
        print*,'avg z_o',sum(z_o(:,1))/dble(n_ll) + 1.d0-rlenz/2.d0
        print*,'delta z_o',maxval(z_o_c)-minval(z_o_c)
      endif
c        if(itime.eq.9) stop
c
      return
      END SUBROUTINE update_h
cc
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
