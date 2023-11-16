      SUBROUTINE comp_gradient(val,dval)
#include "common.inc"

      INTEGER l
      REAL*8 val(n_ll,num_p),dval(n_ll,num_p)

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        dval(l,1)=(val(l,2) - val(l,1))/zplus ! Bottom surface
        dval(l,3)=(val(l,4) - val(l,3))/zplus ! Top    surface

        dval(l,5)=(val(l,6) - val(l,5))/rplus ! First  cylinder
        dval(l,7)=(val(l,8) - val(l,7))/rplus ! Second cylinder
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
! val1 is the Nusselt  number for the 4 surfaces
! val2 is the Sherwood number for the 4 surfaces

      INTEGER l,i,j,n,k
      REAL*8 val1(n_ll,num_p), val2(n_ll,num_p)
      REAL*8 T_fro(n_ll,num_p), rho_fro(n_ll,num_p), p_fro(n_ll,num_p)
      REAL*8 k_fro(n_ll,num_p), Latent(n_ll,num_p), massdiff(n_ll,num_p)
      REAL*8 k_air(n_ll,num_p), eps(n_ll,num_p), T_cel(n_ll,num_p)
      REAL*8 rho_air(n_ll,num_p), rho_fro_n(n_ll,num_p)
      REAL*8 T_fro_n(n_ll,num_p)
      REAL*8 const1(n_ll,num_p), const2(n_ll,num_p)
      REAL*8 const3(n_ll,num_p), const4(n_ll,num_p)
      REAL*8 varA(n_ll,num_p), varB(n_ll,num_p)
      REAL*8 varC(n_ll,num_p), varD(n_ll,num_p)
      REAL*8 dSfrhodt(n_ll,num_p), cp(n_ll,num_p)

      REAL*8 ck1, ck2, ck3, ck4 
      REAL*8 ce0, ce1, ce2, ce3, ce4, ce5
      REAL*8 cp0, cp1, cp2, cp3, cp4, cp5, cp6
      REAL*8 p_atm, R_sgc, nu, Ma, Mw, Va, Vw

      REAL*8 temp_2d(nx0,ny0,num_p), d2tdx2(nx0,ny0,num_p)
      REAL*8 d2tdy2(nx0,ny0,num_p)
      REAL*8 d2tdx2_1d(n_ll,num_p) , d2tdy2_1d(n_ll,num_p)
      REAL*8 dQfdt(n_ll,num_p)

c Coefficients to determine air density 
      p_atm = 101.325d3 ! N/m^2
      R_sgc = 287.058d0 ! specific gas constant J/(kg.K) 
      nu    = 1.42d-5   ! kinematic viscosity of air
      Ma    = 28.9d0    ! Molecular weight of dry air
      Mw    = 18.0d0    ! Molecular weight of water vapor
      Va    = 29.9d0    ! Molar volume of dry air
      Vw    = 18.8d0    ! Molar volume of water vapor

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


!      cp        = 2.04d3 ! specific heat of ice
      do n=1,num_p,2
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c Variables evaluated at old time step
        T_fro(l,n)   = T_wall + (T_inf-T_wall)*tt_p(l,n)!Frost surface temperature
c Frost surface temperature cannot exceed melting temperature
        if (T_fro(l,n).gt.273.15) then
          T_fro(l,n) = 273.15 
        endif

c        T_fro(l,n) = max(T_wall,T_fro(l,n)) 
c Temperature cannot go below wall temperature (Case when Nu<0)
        rho_air(l,n) = p_atm/(R_sgc*T_fro(l,n)) ! rho=p/(R*T)
        k_air(l,n)   = ck1*T_fro(l,n)**3 + ck2*T_fro(l,n)**2
     &                +ck3*T_fro(l,n)    + ck4
        rho_fro(l,n) = 650.d0*dexp(0.227d0*(T_fro(l,n)-273.15d0)) !Mago & Sherif 
        k_fro(l,n)   = 2.4248d-2 + 7.2311d-4*rho_fro(l,n) 
     &                         + 1.1830d-6*rho_fro(l,n)**2
        Latent(l,n)  = -0.1083d0*(1.8d0*(T_fro(l,n)-273.15d0)+32.d0)+2833.d0
        massdiff(l,n)= 0.23d-4 ! We could find a fit for this as well
c        massdiff(l,n)= 4.357d-2*T_fro(l,n)**1.5d0*
c     &        ((dsqrt(1.d0/Ma+1.d0/Mw)
c     &            /(p_atm*(Va**(1.d0/3.d0)+Vw**(1.d0/3.d0))**2)))
        cp(l,n) = (1.85d5 + 6.89d3 * T_fro(l,n)) / rho_fro(l,n)! specific heat of ice



c        tt_p(l,1) = ((z_o(l,1)+z_o_avg-rlenz/2.d0+1.d0) 
c     &               * (k_air(l)/k_fro(l)) * val1(l) + 
c     &               ((z_o(l,1)+z_o_avg-rlenz/2.d0+1.d0) * Y_inf 
c     &               * rho_air(l) * Latent(l) * massdiff(l) / (k_fro(l)
c     &               * (T_inf-T_wall))) * val2(l)) 

      enddo ! l=1,n_ll
C$OMP END PARALLEL DO
      enddo ! n=1,num_p,2

      if (irkk.eq.3) then
c 1 index to 2 index representation
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
          do i=1,nx
            temp_2d(i,j,1) = tt_p((i-1)*ny+j,1)
            temp_2d(i,j,3) = tt_p((i-1)*ny+j,3)
          enddo
        enddo
C$OMP END PARALLEL DO

        call d2dx2fd(temp_2d,d2tdx2)
        call d2dy2fd(temp_2d,d2tdy2)

        n=1 ! Advance temperature for bottom surfaces
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,n_ll
c          dSfrhodt(l,n) = ((z_o(l,1)+z_o_avg-rlenz/2.d0+1.d0)*rho_fro(l,n)-
c     &                    z_o_nm1(l,n)*rho_fro_nm1(l,n))/dt 
          varA(l,n) = cp(l,n)*(z_o(l,n)+1.d0)*rho_fro(l,n) * re*nu/2.d0
          varB(l,n) = k_air(l,n)*val1(l,n)  
     &                + val2(l,n)*Y_inf*rho_air(l,n)*Latent(l,n)
     &                  *massdiff(l,n)/(T_inf-T_wall)
!     &                - cp*re*nu*T_wall*dSfrhodt(l,n)/(T_inf-T_wall,n) 
     &                
          varC(l,n) = k_fro(l,n)/(z_o(l,n)+1.d0)
!     &                +0.5d0*cp*re*nu*dSfrhodt(l,n)

          varD(l,n) = 0.5d0*k_fro(l,n)*(z_o(l,n)+1.d0)

          tt_p(l,n) = (varA(l,n)*tt_p(l,1) + varB(l,n)*AF*dt 
     &                                     + varD(l,n)*AF*dt
     &                *(d2tdx2(1+int((l-1)/ny),1+mod(l-1,ny),n) 
     &                 +d2tdy2(1+int((l-1)/ny),1+mod(l-1,ny),n)))
     &               /(varA(l,n)           + varC(l,n)*AF*dt)
        enddo
C$OMP END PARALLEL DO

        n=3 ! Advance temperature for top    surfaces
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,n_ll
c          dSfrhodt(l,n) = ((z_o(l,1)+z_o_avg-rlenz/2.d0+1.d0)*rho_fro(l,n)-
c     &                    z_o_nm1(l,n)*rho_fro_nm1(l,n))/dt 
          varA(l,n) = cp(l,n)*(-z_o(l,n)+1.d0)*rho_fro(l,n) * re*nu/2.d0
          varB(l,n) = k_air(l,n)*val1(l,n)  
     &                + val2(l,n)*Y_inf*rho_air(l,n)*Latent(l,n)
     &                  *massdiff(l,n)/(T_inf-T_wall)
!     &                - cp*re*nu*T_wall*dSfrhodt(l,n)/(T_inf-T_wall,n) 
     &                
          varC(l,n) = k_fro(l,n)/(-z_o(l,n)+1.d0)
!     &                +0.5d0*cp*re*nu*dSfrhodt(l,n)

          varD(l,n) = 0.5d0*k_fro(l,n)*(-z_o(l,n)+1.d0)

          tt_p(l,n) = (varA(l,n)*tt_p(l,n) + varB(l,n)*AF*dt 
     &                                     + varD(l,n)*AF*dt
     &                *(d2tdx2(1+int((l-1)/ny),1+mod(l-1,ny),n) 
     &                 +d2tdy2(1+int((l-1)/ny),1+mod(l-1,ny),n)))
     &               /(varA(l,n)           + varC(l,n)*AF*dt)
        enddo
C$OMP END PARALLEL DO
! Advance temperature for cylinders
        n=5
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,n_ll
          varA(l,n) = cp(l,n)*(Rpillar+0.5d0*rf(l,1))*rf(l,1)
     &                *rho_fro(l,n)*0.5d0*re*nu*(T_inf-T_wall)
          
          varB(l,n) = (Rpillar+rf(l,1))*(T_inf-T_wall)*k_air(l,n)
     &                 *val1(l,n)
     &                + (Y_inf-Y_wall)*rho_air(l,n)*Latent(l,n)
     &                  *massdiff(l,n)*(Rpillar+rf(l,1))
     &                 *val2(l,n)

          varC(l,n) = k_fro(l,n)*(Rpillar+rf(l,1))*(T_inf-T_wall)
     &                 /rf(l,1)

          tt_p(l,n) = (tt_p(l,n)+AF*dt*varB(l,n)/varA(l,n))
     &               /(1.d0+AF*dt*varC(l,n)/varA(l,n))
        enddo
C$OMP END PARALLEL DO
        n=7
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,n_ll
          varA(l,n) = cp(l,n)*(Rpillar+0.5d0*rf(l,1))*rf(l,1)
     &                *rho_fro(l,n)*0.5d0*re*nu*(T_inf-T_wall)

          varB(l,n) = (Rpillar+rf(l,1))*(T_inf-T_wall)*k_air(l,n)
     &                 *val1(l,n)
     &                + (Y_inf-Y_wall)*rho_air(l,n)*Latent(l,n)
     &                  *massdiff(l,n)*(Rpillar+rf(l,1))
     &                 *val2(l,n)

          varC(l,n) = k_fro(l,n)*(Rpillar+rf(l,1))*(T_inf-T_wall)
     &                 /rf(l,1)

          tt_p(l,n) = (tt_p(l,n)+AF*dt*varB(l,n)/varA(l,n))
     &               /(1.d0+AF*dt*varC(l,n)/varA(l,n))
        enddo
C$OMP END PARALLEL DO
      endif !if (irkk.eq.3) then
c
      do n=1,num_p,2
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c Compute temperature and density at new time step
        T_fro_n(l,n)   = T_wall + (T_inf-T_wall)*tt_p(l,n)!Air temp at frost layer 
c        T_fro_n(l) = max(T_wall,T_fro_n(l)) 
c Temperature cannot go below wall temperature (Case when Nu<0)
        rho_fro_n(l,n) = 650.d0*dexp(0.227d0*(T_fro_n(l,n)-273.15d0)) !Mago & Sherif 
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
        p_fro(l,n) = dexp(cp0/T_fro(l,n)    +cp1 +cp2*T_fro(l,n) 
     &                   +cp3*T_fro(l,n)**2
     &                   +cp4*T_fro(l,n)**3 +cp5*T_fro(l,n)**4 
     &                 +cp6*dlog(T_fro(l,n))) 

        T_cel(l,n) = T_fro(l,n) - 273.15d0
        eps(l,n) = ce0 +ce1*T_cel(l,n)    + ce2*T_cel(l,n)**2 
     &                 +ce3*T_cel(l,n)**3
     &                 +ce4*T_cel(l,n)**4 + ce5*T_cel(l,n)**5 

        tt2_p(l,n) = eps(l,n) * (p_fro(l,n)/p_atm) / Y_inf
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      enddo ! l=1,n_ll
C$OMP END PARALLEL DO
      enddo ! n=1,num_p,2

        print*,itime,irkk,tt_p(1,1)
        print*,itime,irkk,tt_p(1,3)
        print*,itime,irkk,tt_p(1,5)
        print*,itime,irkk,tt_p(1,7)
        print*,itime,irkk,tt2_p(1,1)
        print*,itime,irkk,tt2_p(1,3)
        print*,itime,irkk,tt2_p(1,5)
        print*,itime,irkk,tt2_p(1,7)
c
c      if(mod(itime,oi_mean).eq.0.and.irkk.eq.3) then
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
c        do l=1,n_ll
c          dQfdt(l) = (varB(l) - varC(l)*tt_p(l,1) + varD(l)*
c     &                 (d2tdx2_1d(l) + d2tdy2_1d(l))) / varA(l)
c
cc          write(26,*)l,dQfdt(l)
c        enddo
cC$OMP END PARALLEL DO
cc
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
c        do j=1,ny
c          do i=1,nx
c            d2tdx2_1d((i-1)*ny+j) = d2tdx2(i,j)
c            d2tdy2_1d((i-1)*ny+j) = d2tdy2(i,j)
c          enddo
c        enddo
cC$OMP END PARALLEL DO
c
c! Vapor fraction 
c         write(137,1110) itime, sum(tt2_p(:,1)   ) / dble(n_ll)
c     &                      ,sum(tt2_p(:,1)**2) / dble(n_ll)
c     &                      ,maxval(tt2_p(:,1))
c     &                      ,minval(tt2_p(:,1))
c! Temperature 
c         write(117,1110) itime, sum(tt_p(:,1)   ) / dble(n_ll)
c     &                      ,sum(tt_p(:,1)**2) / dble(n_ll)
c     &                      ,maxval(tt_p(:,1))
c     &                      ,minval(tt_p(:,1))
c! Nusselt
c         write(118,1111) itime, sum(val1   ) / dble(n_ll)
c     &                      ,sum(val1**2) / dble(n_ll)
c     &                      ,maxval(val1)
c     &                      ,minval(val1)
c     &                      ,sum(val1/(1.d0-tt_p(:,1))) / dble(n_ll)
c! Sherwood
c         write(119,1111) itime, sum(val2   ) / dble(n_ll)
c     &                      ,sum(val2**2) / dble(n_ll)
c     &                      ,maxval(val2)
c     &                      ,minval(val2) 
c     &                      ,sum(val2/(1.d0-tt2_p(:,1))) / dble(n_ll)
c! Density  
c         write(120,1110) itime, sum(rho_fro_n   ) / dble(n_ll)
c     &                      ,sum(rho_fro_n**2) / dble(n_ll)
c     &                      ,maxval(rho_fro_n)
c     &                      ,minval(rho_fro_n) 
c! dQfdt_A  
c         write(121,1110) itime, sum(varA   ) / dble(n_ll)
c     &                      ,sum(varA**2) / dble(n_ll)
c     &                      ,maxval(varA)
c     &                      ,minval(varA) 
c! dQfdt_B  
c         write(122,1110) itime, sum(varB   ) / dble(n_ll)
c     &                      ,sum(varB**2) / dble(n_ll)
c     &                      ,maxval(varB)
c     &                      ,minval(varB) 
c! dQfdt_C
c         write(125,1110) itime, sum(varC   ) / dble(n_ll)
c     &                      ,sum(varC**2) / dble(n_ll)
c     &                      ,maxval(varC)
c     &                      ,minval(varC) 
c! dQfdt_C*Theta_f
c         write(126,1110) itime, sum(varC*tt_p(:,1)     ) / dble(n_ll)
c     &                      ,sum((varC*tt_p(:,1))**2) / dble(n_ll)
c     &                      ,maxval(varC*tt_p(:,1))
c     &                      ,minval(varC*tt_p(:,1))
c! dQfdt_D
c         write(127,1110) itime, sum(varD   ) / dble(n_ll)
c     &                      ,sum(varD**2) / dble(n_ll)
c     &                      ,maxval(varD)
c     &                      ,minval(varD) 
c! dQfdt_DE
c         write(128,1110)itime,sum(varD*(d2tdx2_1d+d2tdy2_1d))/dble(n_ll)
c     &                  ,sum((varD*(d2tdx2_1d+d2tdy2_1d))**2)/dble(n_ll)
c     &                      ,maxval(varD*(d2tdx2_1d+d2tdy2_1d))
c     &                      ,minval(varD*(d2tdx2_1d+d2tdy2_1d)) 
c! dQfdt_E1 
c         write(129,1110) itime, sum(dabs(d2tdx2_1d)) / dble(n_ll)
c     &                      ,sum(d2tdx2_1d**2) / dble(n_ll)
c     &                      ,maxval(d2tdx2_1d)
c     &                      ,minval(d2tdx2_1d) 
c! dQfdt_E2
c         write(130,1110) itime, sum(dabs(d2tdy2_1d)) / dble(n_ll)
c     &                      ,sum(d2tdy2_1d**2) / dble(n_ll)
c     &                      ,maxval(d2tdy2_1d)
c     &                      ,minval(d2tdy2_1d) 
c! dQfdt
c         write(131,1110) itime, sum(dQfdt   ) / dble(n_ll)
c     &                      ,sum(dQfdt**2) / dble(n_ll)
c     &                      ,maxval(dQfdt)
c     &                      ,minval(dQfdt)
c      endif

      call update_h(rho_air,rho_fro,rho_fro_n,val2)
      call lagrangian_pts_gen
c
1110  format(I8,4E16.8)
1111  format(I8,5E16.8)
c
      return
      END SUBROUTINE fro_boundary
cc
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
      subroutine update_h(rho_air,rho_fro,rho_fro_n,val2)
#include "common.inc"
      REAL*8 exrc(3), exzc(3) 
      REAL*8 rho_air(n_ll,num_p), rho_fro(n_ll,num_p), val2(n_ll,num_p)
      REAL*8 rho_fro_n(n_ll,num_p)
      REAL*8 z_o_old(n_ll),rf_old(n_ll)
      REAL*8 z_o_2d(nx0,ny0), d2zdx2(nx0,ny0), d2zdy2(nx0,ny0)
      REAL*8 rf_2d(nx0,ny0), d2rdq2(nxyp,nzp),d2rdz2(nxyp,nzp)
      INTEGER l,i,j,k,n
      REAL*8 d2zdx2_1d(n_ll) , d2zdy2_1d(n_ll), dSfdt(n_ll)
      REAL*8 dSfdt1(n_ll), dSfdt2(n_ll)
      REAL*8 dSfdt1_avg

c      diff_coef = 1.d-7

      exrc(1) =  0.d0
      exrc(2) = -5.d0/9.d0
      exrc(3) = -153.d0/128.d0

      exzc(1) = 1.d0/3.d0
      exzc(2) = 15.d0/16.d0
      exzc(3) = 8.d0/15.d0

      n=1 ! Bottom surface
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        z_o_old(l) = z_o(l,n)
      enddo
C$OMP END PARALLEL DO

c 1 index to 2 indices representation
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j=1,ny
        do i=1,nx
          z_o_2d(i,j) = z_o((i-1)*ny+j,n)
        enddo
      enddo
C$OMP END PARALLEL DO

      call d2dx2fd0(z_o_2d,d2zdx2)
      call d2dy2fd0(z_o_2d,d2zdy2)

      dSfdt1_avg = sum((rho_air(:,n)/rho_fro(:,n))
     &                * (Y_inf/(re*pr2)) * val2(:,n))
     &                / dble(nx*ny)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        exner_rhs(l,n) = exrc(irkk)*exner_rhs(l,n)
     &   +AF*dt*((rho_air(l,n)/rho_fro(l,n))*(Y_inf/(re*pr2))*val2(l,n)
     &   +diff_coef*dSfdt1_avg*d2zdx2(1+int((l-1)/ny),1+mod(l-1,ny))
     &   +diff_coef*dSfdt1_avg*d2zdy2(1+int((l-1)/ny),1+mod(l-1,ny)))

        z_o(l,n) =  z_o(l,n) + exzc(irkk)*exner_rhs(l,n)
      enddo
C$OMP END PARALLEL DO

c Compaction
      if(irkk.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,n_ll
          z_o(l,n) = z_o(l,n)-(1.d0-rho_fro(l,n)/rho_fro_n(l,n))
     &           *(z_o_old(l) + 1.d0)
        enddo
C$OMP END PARALLEL DO
      endif

c      
      n=3 ! Bottom surface
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        z_o_old(l) = z_o(l,n) 
      enddo
C$OMP END PARALLEL DO

c 1 index to 2 indices representation
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j=1,ny
        do i=1,nx
          z_o_2d(i,j) = z_o((i-1)*ny+j,n)
        enddo
      enddo
C$OMP END PARALLEL DO
    
      call d2dx2fd0(z_o_2d,d2zdx2)
      call d2dy2fd0(z_o_2d,d2zdy2)

      dSfdt1_avg = sum((rho_air(:,n)/rho_fro(:,n)) 
     &                * (Y_inf/(re*pr2)) * val2(:,n))
     &                / dble(nx*ny)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,n_ll
        exner_rhs(l,n) = exrc(irkk)*exner_rhs(l,n) 
     &   +AF*dt*((rho_air(l,n)/rho_fro(l,n))*(Y_inf/(re*pr2))*val2(l,n)
     &   +diff_coef*dSfdt1_avg*d2zdx2(1+int((l-1)/ny),1+mod(l-1,ny))
     &   +diff_coef*dSfdt1_avg*d2zdy2(1+int((l-1)/ny),1+mod(l-1,ny)))

        z_o(l,n) =  z_o(l,n) - exzc(irkk)*exner_rhs(l,n)
      enddo
C$OMP END PARALLEL DO

c Compaction
      if(irkk.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,n_ll
          z_o(l,n) = z_o(l,n)-(1.d0-rho_fro(l,n)/rho_fro_n(l,n))
     &           *(1.d0 - z_o_old(l))
        enddo
C$OMP END PARALLEL DO
      endif
c
      n=5 ! one of the cylinders
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,nxyp*nzp
        rf_old(l) = rf(l,1)
      enddo
C$OMP END PARALLEL DO

c 1 index to 2 indices representation
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j=1,nzp
        do i=1,nxyp
          rf_2d(i,j) = rf((i-1)*nzp+j,1)
        enddo
      enddo
C$OMP END PARALLEL DO

      call d2dq2fd0(rf_2d,d2rdq2)
      call d2dz2fd0(rf_2d,d2rdz2)

      dSfdt1_avg = sum((rho_air(:,n)/rho_fro(:,n))
     &                * (Y_inf/(re*pr2)) * val2(:,n))
     &                / dble(nxyp*nzp)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,nxyp*nzp
        exner_rhs(l,n) = exrc(irkk)*exner_rhs(l,n)
     &   +AF*dt*((rho_air(l,n)/rho_fro(l,n))*(Y_inf/(re*pr2))*val2(l,n)
     &   +diff_coef*dSfdt1_avg*d2rdq2(1+int((l-1)/nzp),1+mod(l-1,nzp))
     &   +diff_coef*dSfdt1_avg*d2rdz2(1+int((l-1)/nzp),1+mod(l-1,nzp)))

        rf(l,1) =  rf(l,1) + exzc(irkk)*exner_rhs(l,n)
      enddo
C$OMP END PARALLEL DO

c Compaction
      if(irkk.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,nxyp*nzp
          rf(l,1) = rf(l,1)-(1.d0-rho_fro(l,n)/rho_fro_n(l,n))
     &           *rf_old(l)
        enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l,k)
      do k=1,nzp
        do l=1,nxyp
         x_o((l-1)*nzp+k,n)=(Rpillar+rf(l,1))*dcos(dble(l-1)*deltaq)+x01
          if(x_o( (l-1)*nzp+k,n).lt.0) then
                   x_o( (l-1)*nzp+k,n) = x_o( (l-1)*nzp+k,n) + rlenx
           elseif(x_o( (l-1)*nzp+k,n).gt.rlenx) then
                   x_o( (l-1)*nzp+k,n) = x_o( (l-1)*nzp+k,n) - rlenx
           endif

         y_o((l-1)*nzp+k,n)=(Rpillar+rf(l,1))*dsin(dble(l-1)*deltaq)+y01
          if(y_o( (l-1)*nzp+k,n).lt.0) then
                  y_o( (l-1)*nzp+k,n) = y_o( (l-1)*nzp+k,n) + rleny
          elseif(y_o( (l-1)*nzp+k,n).gt.rleny) then
                  y_o( (l-1)*nzp+k,n) = y_o( (l-1)*nzp+k,n) - rleny
          endif
        enddo
      enddo
C$OMP END PARALLEL DO

      endif !irkk.eq.3
c
      n=7 ! one of the cylinders 
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,nxyp*nzp
        rf_old(l) = rf(l,2)
      enddo
C$OMP END PARALLEL DO

c 1 index to 2 indices representation
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j=1,nzp 
        do i=1,nxyp
          rf_2d(i,j) = rf((i-1)*nzp+j,2)
        enddo
      enddo
C$OMP END PARALLEL DO

      call d2dq2fd0(rf_2d,d2rdq2)
      call d2dz2fd0(rf_2d,d2rdz2)

      dSfdt1_avg = sum((rho_air(:,n)/rho_fro(:,n))
     &                * (Y_inf/(re*pr2)) * val2(:,n))
     &                / dble(nxyp*nzp)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
      do l=1,nxyp*nzp
        exner_rhs(l,n) = exrc(irkk)*exner_rhs(l,n)
     &   +AF*dt*((rho_air(l,n)/rho_fro(l,n))*(Y_inf/(re*pr2))*val2(l,n)
     &   +diff_coef*dSfdt1_avg*d2rdq2(1+int((l-1)/nzp),1+mod(l-1,nzp))
     &   +diff_coef*dSfdt1_avg*d2rdz2(1+int((l-1)/nzp),1+mod(l-1,nzp)))

        rf(l,2) =  rf(l,2) + exzc(irkk)*exner_rhs(l,n)
      enddo
C$OMP END PARALLEL DO

c Compaction
      if(irkk.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
        do l=1,nxyp*nzp
          rf(l,2) = rf(l,2)-(1.d0-rho_fro(l,n)/rho_fro_n(l,n))
     &           *rf_old(l)
        enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l,k)
      do k=1,nzp
        do l=1,nxyp
         x_o((l-1)*nzp+k,n)=(Rpillar+rf(l,2))*dcos(dble(l-1)*deltaq)+x02
          if(x_o( (l-1)*nzp+k,n).lt.0) then
                   x_o( (l-1)*nzp+k,n) = x_o( (l-1)*nzp+k,n) + rlenx
           elseif(x_o( (l-1)*nzp+k,n).gt.rlenx) then
                   x_o( (l-1)*nzp+k,n) = x_o( (l-1)*nzp+k,n) - rlenx
           endif

         y_o((l-1)*nzp+k,n)=(Rpillar+rf(l,2))*dsin(dble(l-1)*deltaq)+y02
          if(y_o( (l-1)*nzp+k,n).lt.0) then
                  y_o( (l-1)*nzp+k,n) = y_o( (l-1)*nzp+k,n) + rleny
          elseif(y_o( (l-1)*nzp+k,n).gt.rleny) then
                  y_o( (l-1)*nzp+k,n) = y_o( (l-1)*nzp+k,n) - rleny
          endif

        enddo
      enddo
C$OMP END PARALLEL DO
 
      endif !irkk.eq.3

!      if(mod(itime,oi_mean).eq.0.and.irkk.eq.3) then
!c
!C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
!        do j=1,ny
!          do i=1,nx
!            d2zdx2_1d((i-1)*ny+j) = d2zdx2(i,j)
!            d2zdy2_1d((i-1)*ny+j) = d2zdy2(i,j)
!          enddo
!        enddo
!C$OMP END PARALLEL DO
!
!C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(l)
!        do l=1,n_ll
!          dSfdt1(l) = (rho_air(l)/rho_fro(l))*(Y_inf/(re*pr2))*val2(l)
!          dSfdt2(l) = diff_coef*dsfdt1_avg*(d2zdx2_1d(l)+d2zdy2_1d(l))
!          dSfdt(l)  = dSfdt1(l) + dSfdt2(l)
!        enddo
!C$OMP END PARALLEL DO
!
!! Sf
!         write(116,1110) itime, z_o_avg
!     &    ,sum((z_o_c+1.d0-rlenz/2.d0)**2) / dble(n_ll)
!     &    ,maxval(z_o_c)+1.d0-rlenz/2.d0
!     &    ,minval(z_o_c)+1.d0-rlenz/2.d0
!! dSfdt
!         write(132,1110) itime, sum(dSfdt   ) / dble(n_ll)
!     &                      ,sum(dSfdt**2)    / dble(n_ll)
!     &                      ,maxval(dSfdt)
!     &                      ,minval(dSfdt)
!! dSfdt1
!         write(133,1110) itime, sum(dSfdt1   ) / dble(n_ll)
!     &                      ,sum(dSfdt1**2)    / dble(n_ll)
!     &                      ,maxval(dSfdt1)
!     &                      ,minval(dSfdt1)
!! dSfdt2
!         write(134,1110) itime, sum(dabs(dSfdt2))/ dble(n_ll)
!     &                      ,sum(dSfdt2**2)      / dble(n_ll)
!     &                      ,maxval(dSfdt2)
!     &                      ,minval(dSfdt2)
!! d2Sdx2 same as d2zdx2
!         write(135,1110) itime, sum(dabs(d2zdx2_1d)) / dble(n_ll)
!     &                      ,sum(d2zdx2_1d**2)    / dble(n_ll)
!     &                      ,maxval(d2zdx2_1d)
!     &                      ,minval(d2zdx2_1d)
!! d2Sdy2 same as d2zdy2
!         write(136,1110) itime, sum(dabs(d2zdy2_1d)) / dble(n_ll)
!     &                      ,sum(d2zdy2_1d**2)    / dble(n_ll)
!     &                      ,maxval(d2zdy2_1d)
!     &                      ,minval(d2zdy2_1d)
!c
!      endif
c
1110  format(I8,4E16.8)
c
      return
      END SUBROUTINE update_h
cc
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
