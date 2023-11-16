      SUBROUTINE ibm(cfnl2,cfdif)
#include "common.inc"
c
! This subroutine creates the forcing function (forcing_xyz)
! that will be added to the rhs of the helmholtz equation 
! to account for the immersed body (based on H. Lee)
!
      INTEGER n,l,i,j,k,i1,i2,j1,j2,k1,k2
      REAL*8 u_tilde(nx0,ny0,nz),v_tilde(nx0,ny0,nz)
      REAL*8 w_tilde(nx0,ny0,nz)
      REAL*8 up_tilde(n_ll,num_p),vp_tilde(n_ll,num_p)
      REAL*8 wp_tilde(n_ll,num_p)
      REAL*8 u_forcing(n_ll,num_p),v_forcing(n_ll,num_p)
      REAL*8 w_forcing(n_ll,num_p)
      REAL*8 cfnl2,cfdif
      REAL*8 p_dum(nz,ny,nx),dp_dum(nz,ny,nx)

      INTEGER i_ddf,j_ddf,i_t,j_t,i_f,j_f

      call begin_timer(67)
c--------------------------------------------------------------------
c initializing the variables to zero upon every entry for 
c a new time step
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l,n)
      do n=1,num_p,2
        do l=1,n_l(n)
          up_tilde(l,n)=0.d0
          vp_tilde(l,n)=0.d0
          wp_tilde(l,n)=0.d0
        enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
         do j=1,ny0
            do i=1,nx0
               forcing_x(i,j,k)=0.d0
               forcing_y(i,j,k)=0.d0
               forcing_z(i,j,k)=0.d0
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO

      call ddx_f_new(p,dpx)
      call ddy_f_new(p,dpy)
c      call ddz(p,dpz)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do i=1,nx
         do k=1,nz
            do j=1,ny
               p_dum(k,j,i)=p(i,j,k)
            enddo
         enddo
         call dgemm('N','N',nz,ny,nz,1.d0,chbd1,nz0,
     &          p_dum(1,1,i),nz,0.d0,dp_dum(1,1,i),nz)
         do k=1,nz
            do j=1,ny
               dpz(i,j,k)=dp_dum(k,j,i)
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO


! u_tilde is the actual velocity 
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
         do j=1,ny0
            do i=1,nx0
               u_tilde(i,j,k)=u(i,j,k)+
     &           cfnl2*h1(i,j,k)+2.d0*cfdif*ru(i,j,k)
     &                         -2.d0*cfdif*dpx(i,j,k)
               v_tilde(i,j,k)=v(i,j,k)+
     &           cfnl2*h2(i,j,k)+2.d0*cfdif*rv(i,j,k)
     &                         -2.d0*cfdif*dpy(i,j,k)
               w_tilde(i,j,k)=w(i,j,k)+
     &           cfnl2*h3(i,j,k)+2.d0*cfdif*rw(i,j,k)
     &                         -2.d0*cfdif*dpz(i,j,k)
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
! This section follows the flow solver of M. Uhlmann (2005) to
! define the forcing function.
!
! Transforming the variables to real space
      call fft_f2r_2d_new(u_tilde)
      call fft_f2r_2d_new(v_tilde)
      call fft_f2r_2d_new(w_tilde)
c
      call end_timer(67)
      call begin_timer(68)
cc Transfering u_tilde to the Lagrangian points
C$OMP  PARALLEL DO DEFAULT(SHARED), PRIVATE(n,l,i1,i2,j1,j2,k1,k2
C$OMP& ,i,j,k,i_t,j_t,i_ddf,j_ddf)
      do n=1,num_p,2 ! Skip every other n as velocty gradient not needed
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
                  up_tilde(l,n)=up_tilde(l,n)+u_tilde(i_t,j_t,k)
     &               *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                       (zets(1)-zets(2))/2.d0
                  vp_tilde(l,n)=vp_tilde(l,n)+v_tilde(i_t,j_t,k)
     &                  *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                         (zets(1)-zets(2))/2.d0 
                  wp_tilde(l,n)=wp_tilde(l,n)+w_tilde(i_t,j_t,k) 
     &                  *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                         (zets(1)-zets(2))/2.d0 

                elseif (k .eq. nz) then
                  up_tilde(l,n)=up_tilde(l,n)+u_tilde(i_t,j_t,k)
     &                   *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                     (zets(nz-1)-zets(nz))/2.d0
                  vp_tilde(l,n)=vp_tilde(l,n)+v_tilde(i_t,j_t,k)
     &                   *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                     (zets(nz-1)-zets(nz))/2.d0
                  wp_tilde(l,n)=wp_tilde(l,n)+w_tilde(i_t,j_t,k)
     &                   *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                     (zets(nz-1)-zets(nz))/2.d0

                else
                  up_tilde(l,n)=up_tilde(l,n)+u_tilde(i_t,j_t,k)
     &                   *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                     (zets(k-1)-zets(k+1))/2.d0
                  vp_tilde(l,n)=vp_tilde(l,n)+v_tilde(i_t,j_t,k) 
     &                   *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                     (zets(k-1)-zets(k+1))/2.d0
                  wp_tilde(l,n)=wp_tilde(l,n)+w_tilde(i_t,j_t,k) 
     &                   *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                     (zets(k-1)-zets(k+1))/2.d0

                endif
              enddo
            enddo
          enddo
c
          u_forcing(l,n)=(u_p(l,n)-up_tilde(l,n))!/dt
          v_forcing(l,n)=(v_p(l,n)-vp_tilde(l,n))!/dt
          w_forcing(l,n)=(w_p(l,n)-wp_tilde(l,n))!/dt
        enddo
      enddo
C$OMP END PARALLEL DO

      call end_timer(68)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Caculating the forcing term on the Eulerian grid
      call begin_timer(69)

      do n=1,num_p,2
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
              forcing_x(i_f,j_f,k)=forcing_x(i_f,j_f,k)+u_forcing(l,n)
     &                      *ddf(l,n,i_ddf,j_ddf,k-k1+1)*dv_l(l,n)
              forcing_y(i_f,j_f,k)=forcing_y(i_f,j_f,k)+v_forcing(l,n)
     &                      *ddf(l,n,i_ddf,j_ddf,k-k1+1)*dv_l(l,n)
              forcing_z(i_f,j_f,k)=forcing_z(i_f,j_f,k)+w_forcing(l,n)
     &                      *ddf(l,n,i_ddf,j_ddf,k-k1+1)*dv_l(l,n)
              enddo
            enddo
          enddo
        enddo
cC$OMP END PARALLEL DO
      enddo

      call end_timer(69)
100   format(3(1x,I8),6(3x,g15.8))
c     
c      if( mod((itime-istart),oi_chan/1).eq.0) then
        if(irkk .eq. 1) then
          drag=0.d0
          lift=0.d0
        endif
        do l=1,n_l(1)
           drag=drag+v_forcing(l,1)*dv_l(l,1)
           lift=lift+w_forcing(l,1)*dv_l(l,1)
        enddo
        if(irkk .eq. 3) then
           write(115,*) itime, drag, lift
        endif
c      endif


cc Transforming the forcing function to fourier space
      call fft_r2f_2d_new(forcing_x)
      call fft_r2f_2d_new(forcing_y)
      call fft_r2f_2d_new(forcing_z)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      END SUBROUTINE ibm
