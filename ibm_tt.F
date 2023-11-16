      SUBROUTINE ibm_tt(cfnl2,cfdif)
#include "common.inc"
c
c This subroutine creates the forcing term for the concentration
c with neumann boundary conditions at the immersed boundary
c
      INTEGER n,l,i,j,k,i1,i2,j1,j2,k1,k2
      REAL*8  cfnl2,cfdif
      REAL*8  dzttp_tilde(n_l_max,num_p),dyttp_tilde(n_l_max,num_p)
      REAL*8  tt_forcing(n_l_max,num_p)
      REAL*8  tt_tilde(nx0,ny0,nz)
c      REAL*8  htt_dum(nx0,ny0,nz),rtt_dum(nx0,ny0,nz)
      REAL*8  tt_p(n_l_max,num_p)
      INTEGER i_ddf,j_ddf,i_t,j_t,i_f,j_f
      REAL*8  tt_d_dum(num_p)
cccccccccccccccccccccccccccccccccccccccccccccccccc
      INTEGER lwork,i22,i33,ipiv1(nz)
      PARAMETER (lwork=nz*(nz+4))
      REAL*8  work(lwork)
      REAL*8  tt_tilde_dum(nz,ny,nx),dtt_tilde_dum(nz,ny,nx)
      REAL*8  dztt_tilde(nx0,ny0,nz),dytt_tilde(nx0,ny0,nz)
      REAL*8  dforcing_tt_dum(nz,ny,nx),forcing_tt_dum(nz,ny,nx)
      REAL*8  dforcing_tt(nx,ny,nz)
      REAL*8  tt_c(n_l_max,num_p),dtt_c(n_l_max,num_p)
      REAL*8  ttp_tilde(n_l_max,num_p),wp(n_l_max,num_p)
      REAL*8  dtt_c_dum(n_l_max,num_p)
      REAL*8  dztt_tilde1(nx0,ny0,nz)
c      REAL*8  dforcing_tt_tmp(nx,ny,nz,n_l_max)
cccccccccccccccccccccccccccccccccccccccccccccccccc
      do n=1,num_p
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l)
         do l=1,n_l_max
            dzttp_tilde(l,n)  = 0.d0
            dyttp_tilde(l,n)  = 0.d0
            ttp_tilde  (l,n)  = 0.d0
            wp         (l,n)  = 0.d0
c            tt_forcing (l,n)  = 0.d0
         enddo
C$OMP END PARALLEL DO
      enddo
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
         do j=1,ny0
            do i=1,nx0
               tt_tilde(i,j,k)  = 0.d0
               forcing_tt(i,j,k)= 0.d0
               dforcing_tt(i,j,k)=0.d0
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
c
      if(itime .eq. 1 .and. irkk .eq. 1) then
         do k=1,nz
            do j=1,ny0
               do i=1,nx0
                  forcing_tt(i,j,k)= 0.d0
               enddo
            enddo
         enddo
      endif
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
         do j=1,ny0
            do i=1,nx0
               tt_tilde(i,j,k)=tt(i,j,k)+
     &            cfnl2*htt(i,j,k)+2.d0*cfdif*rtt(i,j,k)
           enddo
         enddo
      enddo
C$OMP END PARALLEL DO
c
      call fft_f2r_2d_new(tt_tilde)
      call fft_f2r_2d_new(w_ibm)
c
      call ddy(tt_tilde,dytt_tilde)
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do i=1,nx
         do k=1,nz
            do j=1,ny
               tt_tilde_dum(k,j,i)=tt_tilde(i,j,k)
            enddo
         enddo
         call dgemm('N','N',nz,ny,nz,1.d0,chbd1,nz0,
     &          tt_tilde_dum(1,1,i),nz,0.d0,dtt_tilde_dum(1,1,i),nz)
         do k=1,nz
            do j=1,ny
               dztt_tilde(i,j,k)=dtt_tilde_dum(k,j,i)
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO
c
      do n=1,num_p
C$OMP  PARALLEL DO DEFAULT(SHARED), PRIVATE(l,i1,i2,j1,j2,k1,k2
C$OMP& ,i,j,k,i_t,i_ddf,j_t,j_ddf)
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
                  if(k .eq. 1) then
                 dzttp_tilde(l,n)=dzttp_tilde(l,n)+dztt_tilde(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                            (zets(1)-zets(2))/2.d0
                 dyttp_tilde(l,n)=dyttp_tilde(l,n)+dytt_tilde(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                            (zets(1)-zets(2))/2.d0
                 ttp_tilde(l,n)=ttp_tilde(l,n)+tt(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                            (zets(1)-zets(2))/2.d0
                 wp(l,n)       =wp(l,n)+w_ibm(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                            (zets(1)-zets(2))/2.d0
                  elseif (k .eq. nz) then
                 dzttp_tilde(l,n)=dzttp_tilde(l,n)+dztt_tilde(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                        (zets(nz-1)-zets(nz))/2.d0
                 dyttp_tilde(l,n)=dyttp_tilde(l,n)+dytt_tilde(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                        (zets(nz-1)-zets(nz))/2.d0
                 ttp_tilde(l,n)=ttp_tilde(l,n)+tt(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                        (zets(nz-1)-zets(nz))/2.d0
                 wp(l,n)       =wp(l,n)+w_ibm(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                        (zets(nz-1)-zets(nz))/2.d0
                  else
                 dzttp_tilde(l,n)=dzttp_tilde(l,n)+dztt_tilde(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                        (zets(k-1)-zets(k+1))/2.d0
                 dyttp_tilde(l,n)=dyttp_tilde(l,n)+dytt_tilde(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                        (zets(k-1)-zets(k+1))/2.d0
                 ttp_tilde(l,n)=ttp_tilde(l,n)+tt(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                        (zets(k-1)-zets(k+1))/2.d0
                 wp(l,n)       =wp(l,n)+w_ibm(i_t,j_t,k)
     &                       *ddf(l,n,i_ddf,j_ddf,k-k1+1)*deltax*deltay*
     &                                        (zets(k-1)-zets(k+1))/2.d0
                  endif
              enddo
             enddo
            enddo
         enddo
C$OMP END PARALLEL DO
      enddo
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l)
      do l=1,n_l(1)
         dtt_c_dum(l,1)= dyttp_tilde(l,1)*dcos(theta_ibm(l,1))
     &                      +dzttp_tilde(l,1)*dsin(theta_ibm(l,1))
         dtt_c(l,1)=(3.d0/2.d0*dtt_c_dum(l,1)-
     &                1.d0/2.d0*dtt_c_dum_n(l,1))*dsin(theta_ibm(l,1))
         dtt_c_dum_n(l,1)=dtt_c_dum(l,1)
c
       tt_c(l,1)=(wp(l,1)-ws)*re*pr*(dcos(theta_ibm(l,1)))**2.d0*
     &         (3.d0/2.d0*ttp_tilde(l,1)-1.d0/2.d0*ttp_tilde_n(l,1))
         ttp_tilde_n(l,1)=ttp_tilde(l,1)
      enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(l)
      do l=1,n_l(1)
         tt_forcing(l,1)=((tt_c(l,1)+dtt_c(l,1))-dzttp_tilde(l,1))/dt
      enddo
C$OMP END PARALLEL DO
c
      do n=1,1
C$OMP  PARALLEL DO DEFAULT(SHARED), PRIVATE(l,i1,i2,j1,j2,k1,k2
C$OMP& ,i,j,k,i_f,j_f,i_ddf,j_ddf)
       do l=1,num_openmp_ibm*openmp_proc
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
               dforcing_tt(i_f,j_f,k)=dforcing_tt(i_f,j_f,k)+
     &            tt_forcing(l,n)*ddf(l,n,i_ddf,j_ddf,k-k1+1)*dv_l(l,n)
            enddo
           enddo
         enddo
       enddo
C$OMP END PARALLEL DO
      enddo
c
      do n=1,1
       do l=num_openmp_ibm*openmp_proc+1,n_l(n)
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
               dforcing_tt(i_f,j_f,k)=dforcing_tt(i_f,j_f,k)+
     &            tt_forcing(l,n)*ddf(l,n,i_ddf,j_ddf,k-k1+1)*dv_l(l,n)
            enddo
           enddo
         enddo
       enddo
      enddo
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      do k=1,nz
c         chbd1_ibm(k_ibm-2,k)=0.d0
c      enddo
c      chbd1_ibm(k_ibm-2,1)=1.d0
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do i=1,nx
        do j=1,ny
          do k=1,nz
            dforcing_tt_dum(k,j,i)=dforcing_tt(i,j,k)
          enddo
c        call dgemm('N','N',5,ny,nz,1.d0,chbd1_ibm(k_ibm-2,1),nz,
c     &                      dforcing_tt_dum(k_ibm-2,1,i),nz,0.d0,
c     &                            forcing_tt_dum(k_ibm-2,1,i),nz)
         call dgemm('N','N',nz,ny,nz,1.d0,chbd1_ibm,nz,
     &                      dforcing_tt_dum(1,1,i),nz,0.d0,
     &                            forcing_tt_dum(1,1,i),nz)
        enddo
      enddo
C$OMP END PARALLEL DO
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do i=1,nx
        do j=1,ny
          do k=k_ibm(i,j)-2,nz
c         do k=1,nz
               forcing_tt(i,j,k)=forcing_tt_dum(k,j,i)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO
cC$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,i,j)
c      do k=k_ibm+5,nz
c         do i=1,nx
c            do j=1,ny
c               forcing_tt(i,j,k)=0.d0
c            enddo
c         enddo
c      enddo
ccC$OMP END PARALLEL DO
c      call fft_f2r_2d_new(u)
c      call fft_f2r_2d_new(v)
c      call fft_f2r_2d_new(w)
c      call fft_f2r_2d_new(tt)
cC$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(k,i,j)
c      do k=k_ibm+2,nz
c         do i=1,nx
c            do j=1,ny
c               tt(i,j,k)=0.d0 !0.5d0*tt(i,j,k_ibm+7)+0.5d0*tt(i,j,k_ibm+7)
cc     &                       *dtanh((zets(k)-zets(k_ibm+2))/delta_l)
c               u(i,j,k)=0.d0 !0.5d0*u(i,j,k_ibm+7)+0.5d0*u(i,j,k_ibm+7)
cc     &                       *dtanh((zets(k)-zets(k_ibm+2))/delta_l)
c               v(i,j,k)=0.d0 !0.5d0*v(i,j,k_ibm+7)+0.5d0*v(i,j,k_ibm+7)
cc     &                       *dtanh((zets(k)-zets(k_ibm+2))/delta_l)
c               w(i,j,k)=0.d0 !0.5d0*w(i,j,k_ibm+7)+0.5d0*w(i,j,k_ibm+7)
cc     &                       *dtanh((zets(k)-zets(k_ibm+2))/delta_l)
c            enddo
c         enddo
c      enddo
cC$OMP END PARALLEL DO
c      call fft_r2f_2d_new(u)
c      call fft_r2f_2d_new(v)
c      call fft_r2f_2d_new(w)
c      call fft_r2f_2d_new(tt)
c
      call fft_r2f_2d_new(forcing_tt)
c
      END SUBROUTINE ibm_tt
