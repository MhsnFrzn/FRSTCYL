      SUBROUTINE ibm_tt_new(cfnl2,cfdif)
#include "common.inc"
c
! This subroutine creates the forcing function (forcing_tt)
! that will be added to the rhs of the helmholtz equation 
! to account for the immersed body (based on H. Lee)
!
      INTEGER n,l,i,j,k,i1,i2,j1,j2,k1,k2
      REAL*8 tt_tilde(nx0,ny0,nz)
      REAL*8 ttp_tilde(n_ll,num_p), tt_forcing(n_ll,num_p)
      REAL*8 cfnl2,cfdif

      INTEGER i_ddf,j_ddf,i_t,j_t,i_f,j_f

      call begin_timer(67)
c--------------------------------------------------------------------
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny0
          do i=1,nx0
            tt_tilde(i,j,k)=tt(i,j,k)+
     &          cfnl2*htt(i,j,k)+2.d0*cfdif*rtt(i,j,k)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
! Taking tt_tilde to real space
      call fft_f2r_2d_new(tt_tilde)
c
      call begin_timer(68)
! We need to loop here over ALL num_p values to interpolate at 2nd layer
      call Euler2Lagrange(1,num_p,1,tt_tilde,ttp_tilde)
      call end_timer(68)
c      if(mod(itime,oi_mean).eq.0.and.irkk.eq.3) then
c        print*,'ttp_tilde(1,1) ',ttp_tilde(1,1)
c        print*,'tt_p(1,1) ',tt_p(1,1)
c        print*,'ttp_tilde(1,2) ',ttp_tilde(1,2)
c        print*,'tt_p(1,2) ',tt_p(1,2)
c      endif
c
      if(FLMOV) call comp_gradient(ttp_tilde,dttp_tildedn)
c      if(mod(itime,oi_mean).eq.0.and.irkk.eq.3) then
c        print*,itime,irkk
c        print*,'tt(1,1,1) = ',tt(1,1,1)
c        print*,'tt2(1,1,1) = ',tt2(1,1,1)
c        print*,'max dttp_tildedn',maxval(dttp_tildedn)
c        print*,'min dttp_tildedn',minval(dttp_tildedn)
c      endif

c
! This section follows the flow solver of M. Uhlmann (2005) to
! define the forcing function.
      do n=1,num_p,2
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l)
         do l=1,n_l(n)
c           tt_p(l,n)   = 0.d0 ! Nadim: This needs to be updated at each stage
           tt_forcing(l,n)=(tt_p(l,n)-ttp_tilde(l,n)) ! Nadim: No longer dividing by dt! 
         enddo
C$OMP END PARALLEL DO
      enddo
c
      call Lagrange2Euler(1,num_p,2,tt_forcing,forcing_tt)
c
      call end_timer(67)
c
100   format(3(1x,I8),6(3x,g15.8))
c     
cc Transforming the forcing function to fourier space
      call fft_r2f_2d_new(forcing_tt)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      END SUBROUTINE ibm_tt_new
c
      SUBROUTINE ibm_tt2_new(cfnl2,cfdif)
#include "common.inc"
c
! This subroutine creates the forcing function (forcing_tt)
! that will be added to the rhs of the helmholtz equation 
! to account for the immersed body (based on H. Lee)
!
      INTEGER n,l,i,j,k,i1,i2,j1,j2,k1,k2
      REAL*8 tt_tilde(nx0,ny0,nz)
      REAL*8 ttp_tilde(n_ll,num_p), tt_forcing(n_ll,num_p)
      REAL*8 cfnl2,cfdif

      INTEGER i_ddf,j_ddf,i_t,j_t,i_f,j_f

      call begin_timer(67)
c--------------------------------------------------------------------
c
! This section follows the flow solver of M. Uhlmann (2005) to
! define the forcing function.
      n=1
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
        do j=1,ny0
          do i=1,nx0
            tt_tilde(i,j,k)=tt2(i,j,k)+
     &          cfnl2*htt2(i,j,k)+2.d0*cfdif*rtt2(i,j,k)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
! Taking tt2_tilde to real space
      call fft_f2r_2d_new(tt_tilde)
c
      call begin_timer(68)
      call Euler2Lagrange(1,num_p,1,tt_tilde,ttp_tilde)
      call end_timer(68)
c      if(mod(itime,oi_mean).eq.0.and.irkk.eq.3) then
c        print*,'tt2p_tilde(1,1) ',ttp_tilde(1,1)
c        print*,'tt2_p(1,1) ',tt2_p(1,1)
c        print*,'tt2p_tilde(1,2) ',ttp_tilde(1,2)
c        print*,'tt2_p(1,2) ',tt2_p(1,2)
c      endif
c  
      if(FLMOV) call comp_gradient(ttp_tilde,dtt2p_tildedn)
c  
      if(FLMOV) call fro_boundary(dttp_tildedn,dtt2p_tildedn)

      do n=1,num_p,2
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l)
         do l=1,n_l(n)
c           if(.not.FLMOV) tt2_p(l,n)   = 0.102324008830544d0  ! Nadim: This needs to be updated at each stage
           tt_forcing(l,n)=(tt2_p(l,n)-ttp_tilde(l,n)) ! Nadim: No longer dividing by dt!
         enddo
C$OMP END PARALLEL DO
      enddo
c  
      call Lagrange2Euler(1,num_p,2,tt_forcing,forcing_tt2)
!
      call end_timer(67)
c
100   format(3(1x,I8),6(3x,g15.8))
c     
cc Transforming the forcing function to fourier space
      call fft_r2f_2d_new(forcing_tt2)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      END SUBROUTINE ibm_tt2_new
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Lagrange2Euler(ni,nf,ninc,Lagrange,Euler)
#include "common.inc"

      INTEGER n,ni,nf,ninc,l,i,j,k,i1,i2,j1,j2,k1,k2
      REAL*8 Euler(nx0,ny0,nz), Lagrange(n_ll,num_p)

      INTEGER i_ddf,j_ddf,i_f,j_f

C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
        do k=1,nz
          do j=1,ny0
            do i=1,nx0
               Euler(i,j,k) = 0.d0
            enddo
          enddo
        enddo
C$OMP END PARALLEL DO

      do n=ni,nf,ninc
cC$OMP  PARALLEL DO DEFAULT(SHARED), PRIVATE(l,i1,i2,j1,j2,k1,k2
cC$OMP& ,i,j,k,i_f,j_f,i_ddf,j_ddf),
cC$OMP&  REDUCTION(+:forcing_x,forcing_y,forcing_z)! Using Reduction is
c about 2 orders of magnitude slower than serial computation for a
c resolution of nxnynz=4x1200x201 n_ll 4840
        do l=1,n_l(n)
          i1=p_iw(l,n)
          i2=p_ie(l,n)
          j1=p_js(l,n)
          j2=p_jn(l,n)
          k1=p_kb(l,n)
          k2=p_kt(l,n)
          do i=i1,i2
            i_ddf=i-i1+1
            if(i .lt. 1) then
              i_f=nx+i
            elseif (i .gt. nx) then
              i_f=i-nx
            else
              i_f=i
            endif
           do j=j1,j2
              j_ddf=j-j1+1
              if(j .lt.1) then
                 j_f=ny+j
              elseif(j .gt. ny) then
                j_f=j-ny
              else
                j_f=j
              endif
              do k=k1,k2
                Euler(i_f,j_f,k) = Euler(i_f,j_f,k) +
     &             Lagrange(l,n)*ddf(l,n,i_ddf,j_ddf,k-k1+1)
     &               *dv_l(l,n)
              enddo
            enddo
          enddo
        enddo
cC$OMP END PARALLEL DO
      enddo
c
      return
      END SUBROUTINE Lagrange2Euler
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE Euler2Lagrange(ni,nf,ninc,Euler,Lagrange)
#include "common.inc"

      INTEGER n,ni,nf,ninc,l,i,j,k,i1,i2,j1,j2,k1,k2
      REAL*8 Euler(nx0,ny0,nz), Lagrange(n_ll,num_p)

      INTEGER i_ddf,j_ddf,i_t,j_t,i_f,j_f

      do n=ni,nf,ninc
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l)
        do l=1,n_l(n)
          Lagrange(l,n)=0.d0
        enddo
C$OMP END PARALLEL DO
        enddo
      do n=ni,nf,ninc
C$OMP  PARALLEL DO DEFAULT(SHARED), PRIVATE(l,i1,i2,j1,j2,k1,k2
C$OMP& ,i,j,k,i_t,j_t,i_ddf,j_ddf)
        do l=1,n_l(n)
          i1=p_iw(l,n)
          i2=p_ie(l,n)
          j1=p_js(l,n)
          j2=p_jn(l,n)
          k1=p_kb(l,n)
          k2=p_kt(l,n)
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
                if (k .eq. 1) then
                  Lagrange(l,n)=Lagrange(l,n)+Euler(i_t,j_t,k)
     &                     *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                          (zets(1)-zets(2))/2.d0
                elseif (k .eq. nz) then
                  Lagrange(l,n)=Lagrange(l,n)+Euler(i_t,j_t,k)
     &                     *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                      (zets(nz-1)-zets(nz))/2.d0

                else
                  Lagrange(l,n)=Lagrange(l,n)+Euler(i_t,j_t,k)
     &                     *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                      (zets(k-1)-zets(k+1))/2.d0
                endif
              enddo
            enddo
          enddo

        enddo
C$OMP END PARALLEL DO
      enddo

      return
      END SUBROUTINE Euler2Lagrange

