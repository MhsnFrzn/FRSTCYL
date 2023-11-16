c     This computes the x- and y- spectra for energy and dissipation.
c     Contributions from the mean profile aren't included.
c     Also outputs the rms energy as a function of z, and the rms
c     average of this over the box
c
c     xrms and yrms should be the same, as should txrms and tyrms.
c     (Both are included as a check.)
c
      SUBROUTINE make_spec_new(uf,espx,dspx,espy,dspy,xrms,txrms,
     &  yrms,tyrms)
#include "common.inc"
      REAL*8 uf(nx0,ny0,nz)
      REAL*8 espx(nxh,nz),dspx(nxh,nz),espy(nyh,nz),dspy(nyh,nz)
      REAL*8 xrms(nz),txrms,yrms(nz),tyrms
c
      INTEGER i,ii,j,jimag,k
      REAL*8 ensum,dssum,temp
c
      if (FLTIMER) call begin_timer(110)
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP&            PRIVATE(i,ii,j,jimag,k,ensum,dssum,temp)
      do k=1,nz
c
c     compute x-spectra
c
        ensum = 0.d0
        dssum = 0.d0
        do j=2,nyh
          jimag=nyh+j
          temp = uf(1,j,k)**2 + uf(1,jimag,k)**2
          ensum = ensum + temp
          dssum = dssum + waveys(j)*temp
        enddo
        espx(1,k) = 2.d0*ensum
        dspx(1,k) = 2.d0*dssum
        xrms(k) = espx(1,k)
        do i=2,nxh
          ii = nx+2-i
          ensum = 0.5d0*(uf(i, 1,k)**2 + uf(i, nyhp,k)**2 +
     &                   uf(ii,1,k)**2 + uf(ii,nyhp,k)**2)
          dssum = 0.d0
          do j=2,nyh
            jimag=nyh+j
            temp = uf(i, j,k)**2 + uf(i, jimag,k)**2 +
     &             uf(ii,j,k)**2 + uf(ii,jimag,k)**2
            ensum = ensum + temp
            dssum = dssum + waveys(j)*temp
          enddo
          espx(i,k) = 2.d0*ensum
          dspx(i,k) = 2.d0*(wavexs(i)*ensum + dssum)
          xrms(k) = xrms(k) + espx(i,k)
        enddo
c
c     compute y-spectra
c
        ensum = 0.d0
        dssum = 0.d0
        do i=2,nxh
          ii = nx+2-i
          temp = uf(i, 1,k)**2 + uf(i, nyhp,k)**2 +
     &           uf(ii,1,k)**2 + uf(ii,nyhp,k)**2
          ensum = ensum + temp
          dssum = dssum + wavexs(i)*temp
        enddo
        espy(1,k) = ensum
        dspy(1,k) = dssum
        yrms(k) = espy(1,k)
        do j=2,nyh
          jimag=nyh+j
          ensum = uf(1,j,k)**2 + uf(1,jimag,k)**2
          dssum = 0.d0
          do i=2,nxh
            ii = nx+2-i
            temp = uf(i, j,k)**2 + uf(i, jimag,k)**2 +
     &             uf(ii,j,k)**2 + uf(ii,jimag,k)**2
            ensum = ensum + temp
            dssum = dssum + wavexs(i)*temp
          enddo
          espy(j,k) = 2.d0*ensum
          dspy(j,k) = 2.d0*(waveys(j)*ensum + dssum)
          yrms(k) = yrms(k) + espy(j,k)
        enddo
      enddo
C$OMP END PARALLEL DO
c
      txrms = 0.d0
      tyrms = 0.d0
      do k=1,nz
        txrms = txrms + xrms(k)*wzet(k)
        tyrms = tyrms + yrms(k)*wzet(k)
        xrms(k) = sqrt(xrms(k))
        yrms(k) = sqrt(yrms(k))
      enddo
      txrms = sqrt(txrms)
      tyrms = sqrt(tyrms)
c
      if (FLTIMER) call end_timer(110)
      return
      END
c
c     This computes Chebyshev energy spectra 
c
      SUBROUTINE make_chb_spec(val,espz)
#include "common.inc"
      REAL*8 val(nx0,ny0,nz)
      REAL*8 espz(nz)
c
      INTEGER i,j,jcomp,k
c
c     First half
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j)
      do j = 1,nyh
c     transpose data
        do k = 1,nz
          call dcopy(nx,val(1,j,k),1,auxt1(k,1,j),nz0)
        enddo
c     compute
        call dgemm('N','N',nz,nx,nz,1.d0,chb,nz0,auxt1(1,1,j),nz0,
     &                                      0.d0,auxt2(1,1,j),nz0)
      enddo
C$OMP END PARALLEL DO
c     
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
       espz(k)=0.d0
       do j=1,nyh 
        do i=1,nx
         espz(k)=espz(k)+auxt2(k,i,j)**2
        enddo
       enddo
       espz(k)=espz(k)/dble(nx*ny)
      enddo
C$OMP END PARALLEL DO
c     Second half
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,jcomp)
      do j = 1,nyh
        jcomp=nyh+j
c     transpose data
        do k = 1,nz
          call dcopy(nx,val(1,jcomp,k),1,auxt1(k,1,j),nz0)
        enddo
c     compute
        call dgemm('N','N',nz,nx,nz,1.d0,chb,nz0,auxt1(1,1,j),nz0,
     &                                      0.d0,auxt2(1,1,j),nz0)
      enddo
C$OMP END PARALLEL DO
c     
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
       do j=1,nyh
        do i=1,nx
         espz(k)=espz(k)+auxt2(k,i,j)**2
        enddo
       enddo
       espz(k)=espz(k)/dble(nx*ny)
      enddo
C$OMP END PARALLEL DO
c
      return
      END 
c
c     This computes the mean as a function of z
c     and the total mean over the box.
c     Uses Fourier space information.
c
      SUBROUTINE make_meanf(ur,mean,tmean)
#include "common.inc"
      REAL*8 ur(nx0,ny0,nz),mean(nz),tmean
c
      INTEGER k
c
      tmean = 0.d0
      do k=1,nz
        mean(k) = ur(1,1,k)
        tmean = tmean + mean(k)*wzet(k)
      enddo
c
      return
      END
c
c     This computes the mean and rms as a function of z
c     and the total mean and rms over the box.
c     Uses physical space information.
c
      SUBROUTINE make_meanp(val,mean,tmean,rms,trms)
#include "common.inc"
      REAL*8 val(nx0,ny0,nz),mean(nz),tmean,rms(nz),trms
c
      INTEGER i,j,k
      REAL*8 s,su
c
      if (FLTIMER) call begin_timer(111)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,s,su)
      do k=1,nz
        su = 0.d0
        do j=1,ny
          s = 0.d0
          do i=1,nx
            s = s + val(i,j,k)
          enddo
          su = su + s
        enddo
        mean(k) = su/dble(nx*ny)
c
        su = 0.d0
        do j=1,ny
          s = 0.d0
          do i=1,nx
            s = s + (val(i,j,k)-mean(k))**2
          enddo
          su = su + s
        enddo
        rms(k) = su/dble(nx*ny)
      enddo
C$OMP END PARALLEL DO
c
      tmean = 0.d0
      trms = 0.d0
      do k=1,nz
        tmean = tmean + mean(k)*wzet(k)
        trms = trms + rms(k)*wzet(k)
        rms(k) = sqrt(rms(k))
      enddo
      trms = sqrt(trms)
c
      if (FLTIMER) call end_timer(111)
      return
      END
c
c     This computes the mean square of a quantity
c
      SUBROUTINE make_msq(val,msq)
#include "common.inc"
      REAL*8 val(nx0,ny0,nz),msq
c
      INTEGER i,j,k
      REAL*8 s,su,msqz(nz)
c
      if (FLTIMER) call begin_timer(112)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,s,su)
      do k=1,nz
        su = 0.d0
        do j=1,ny
          s = 0.d0
          do i=1,nx
            s = s + val(i,j,k)**2
          enddo
          su = su + s
        enddo
        msqz(k) = su/dble(nx*ny)
      enddo
C$OMP END PARALLEL DO
c
      msq = 0.d0
      do k=1,nz
        msq = msq + msqz(k)*wzet(k)
      enddo
c
      if (FLTIMER) call end_timer(112)
      return
      END
c
c     This computes the quantities on the RHS of global balances
c     Uses physical space quantities
c
      SUBROUTINE make_gbalrhs(KE,Ep,diss,ddtEp)
#include "common.inc"
      REAL*8 KE,Ep,diss,ddtEp
c
      INTEGER i,j,k
      REAL*8 mean(nz),tmean,s,su
c
      if (FLTIMER) call begin_timer(115)
c
c     compute KE and 
c
      KE = 0.d0
      call make_msq(u,tmean)
      KE = KE + tmean
      call make_msq(v,tmean)
      KE = KE + tmean
      call make_msq(w,tmean)
      KE = KE + tmean
c 
      KE = KE*rlenx*rleny*rlenz
c
c
c     compute KE dissipation rate
c
      diss = 0.d0
c
      call ddx(u,ru)
      call make_msq(ru,tmean)
      diss = diss + tmean
      call ddy(v,ru)
      call make_msq(ru,tmean)
      diss = diss + tmean
      call ddz(w,ru)
      call make_msq(ru,tmean)
      diss = diss + tmean
c
      call ddy(u,ru)
      call ddx(v,rv)
      call ddz(u,rw)
      call ddx(w,rtt)
      call ddy(w,rp)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       do k=1,nz
        do j=1,ny
         do i=1,nx
           ru(i,j,k)=(ru(i,j,k)+ rv(i,j,k))/2.d0
           rw(i,j,k)=(rw(i,j,k)+rtt(i,j,k))/2.d0
         enddo
        enddo
       enddo
C$OMP END PARALLEL DO
      call ddz(v,rv)
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
       do k=1,nz
        do j=1,ny
         do i=1,nx
            rv(i,j,k)=( rv(i,j,k)+  rp(i,j,k))/2.d0
         enddo
        enddo
       enddo
C$OMP END PARALLEL DO
      call make_msq(ru,tmean)
      diss = diss + 2.d0*tmean
      call make_msq(rw,tmean)
      diss = diss + 2.d0*tmean
      call make_msq(rv,tmean)
      diss = diss + 2.d0*tmean
c
      diss = 2.d0*diss/re*(rlenx*rleny*rlenz)
c
c
c     compute Ep 
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
        mean(k)=0.d0
        do j=1,ny
          do i=1,nx
            mean(k) = mean(k) + (zets(k)+1.d0)*tt(i,j,k)
          enddo
        enddo
        mean(k) = mean(k)*deltax*deltay
      enddo
C$OMP END PARALLEL DO
c
      Ep = 0.d0
      do k=1,nz
        Ep = Ep + mean(k)*wzet(k)*rlenz
      enddo
c
c
c      compute ddtEp
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
      do k=1,nz
        mean(k)=0.d0
        do j=1,ny
          do i=1,nx
            mean(k) = mean(k) + w(i,j,k)*tt(i,j,k)
          enddo
        enddo
        mean(k) = mean(k)*deltax*deltay
      enddo
C$OMP END PARALLEL DO
c
      ddtEp = 0.d0
      do k=1,nz
        ddtEp = ddtEp + mean(k)*wzet(k)*rlenz
      enddo
c
      if(FLSTLE)then
       call ddz(tt,rtt)
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
       do k=1,nz
         mean(k)=0.d0
         do j=1,ny
           do i=1,nx
             mean(k) = mean(k) + (zets(k)+1.d0)*ws*rtt(i,j,k)
           enddo
         enddo
         mean(k) = mean(k)*deltax*deltay
       enddo
C$OMP END PARALLEL DO        
c
       do k=1,nz
         ddtEp = ddtEp + mean(k)*wzet(k)*rlenz
       enddo
      endif
c
      if (FLTIMER) call end_timer(115)
      return
      END
c
c
c     This computes the accumulative means
c     Assumes variable in Fourier space and all reusable variables free
c
      SUBROUTINE make_means
#include "common.inc"
c     
      INTEGER i,j,k,ip,itimeaux
      REAL*8 sum,svm,swm,sttm,su2m,sv2m,sw2m,stt2m,su3m,sv3m,sw3m,stt3m,
     &       su4m,sv4m,sw4m,stt4m,suvm,suwm,suttm,svwm,svttm,swttm,
     &       spm,spwm,su2wm,sv2wm,suuzm,suwxm,svvzm,svwym,swwzm
c
      REAL*8 espx(nxh,nz),dspx(nxh,nz),espy(nyh,nz),dspy(nyh,nz)
      REAL*8 xrms(nz),txrms,yrms(nz),tyrms
c
      REAL*8 ddot
      EXTERNAL ddot
c
      REAL*8 aux,avx,awx,attx,rp_dum(nx0,ny0,nz)
c      INTEGER k3,k5,k10,k20,k50,k100,k200,k300,knz4,k3nz8,knzh
c      k3    = 5
c      k5    = 7  !6
c      k10   = 9  !8
c      k20   = 13 ! 10
c      k50   = 20 !16
c      k100  = 28 !22
c      k200  = 39 !30
c      k300  = 48 !37
c      knz4  = nz/4
c      k3nz8 = 3*nz/8
c      knzh  = nz/2
c
c     produce physical space variables
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k)
      do k=1,nz
         call dcopy(nx0y0,u( 1,1,k),1,ru( 1,1,k),1)
         call dcopy(nx0y0,v( 1,1,k),1,rv( 1,1,k),1)
         call dcopy(nx0y0,w( 1,1,k),1,rw( 1,1,k),1)
         call dcopy(nx0y0,p( 1,1,k),1,rp_dum( 1,1,k),1)
      enddo
C$OMP END PARALLEL DO
      call fft_f2r_2d_new(ru)
      call fft_f2r_2d_new(rv)
      call fft_f2r_2d_new(rw)
      call fft_f2r_2d_new(rp_dum)
c
      if(FLTHRM)then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k)
         do k=1,nz
            call dcopy(nx0y0,tt(1,1,k),1,rtt(1,1,k),1)
         enddo
C$OMP END PARALLEL DO
         call fft_f2r_2d_new(rtt)
      endif
c
c     output time series
cC$OMP PARALLEL DEFAULT(SHARED)
cC$OMP SECTIONS
cC$OMP SECTION
c      write(60,300)rtime,rv(nx/2,1,kts(1)),rv(nx/2,1,kts(2))
c     & ,rv(nx/2,1,kts(3)),rv(nx/2,1,kts(4)),rv(nx/2,1,kts(5))
c     & ,rv(nx/2,1,kts(6)),rv(nx/2,1,kts(7)),rv(nx/2,1,kts(8))
c     & ,rv(nx/2,1,kts(9)),rv(nx/2,1,kts(10)),rv(nx/2,1,kts(11))
c      write(61,300)rtime,ru(nx/2,1,kts(1)),ru(nx/2,1,kts(2))
c     & ,ru(nx/2,1,kts(3)),ru(nx/2,1,kts(4)),ru(nx/2,1,kts(5))
c     & ,ru(nx/2,1,kts(6)),ru(nx/2,1,kts(7)),ru(nx/2,1,kts(8))
c     & ,ru(nx/2,1,kts(9)),ru(nx/2,1,kts(10)),ru(nx/2,1,kts(11))
c      write(62,300)rtime,rw(nx/2,1,kts(1)),rw(nx/2,1,kts(2))
c     & ,rw(nx/2,1,kts(3)),rw(nx/2,1,kts(4)),rw(nx/2,1,kts(5))
c     & ,rw(nx/2,1,kts(6)),rw(nx/2,1,kts(7)),rw(nx/2,1,kts(8))
c     & ,rw(nx/2,1,kts(9)),rw(nx/2,1,kts(10)),rw(nx/2,1,kts(11))
c      write(70,300)rtime,rv(1,1,nz+1-kts(1)),rv(1,1,nz+1-kts(2))
c     & ,rv(1,1,nz+1-kts(3)),rv(1,1,nz+1-kts(4)),rv(1,1,nz+1-kts(5))
c     & ,rv(1,1,nz+1-kts(6)),rv(1,1,nz+1-kts(7)),rv(1,1,nz+1-kts(8))
c     & ,rv(1,1,nz+1-kts(9)),rv(1,1,nz+1-kts(10)),rv(1,1,nz+1-kts(11))
cC$OMP SECTION
c      write(71,300)rtime,ru(1,1,nz+1-kts(1)),ru(1,1,nz+1-kts(2))
c     & ,ru(1,1,nz+1-kts(3)),ru(1,1,nz+1-kts(4)),ru(1,1,nz+1-kts(5))
c     & ,ru(1,1,nz+1-kts(6)),ru(1,1,nz+1-kts(7)),ru(1,1,nz+1-kts(8))
c     & ,ru(1,1,nz+1-kts(9)),ru(1,1,nz+1-kts(10)),ru(1,1,nz+1-kts(11))
c      write(72,300)rtime,rw(1,1,nz+1-kts(1)),rw(1,1,nz+1-kts(2))
c     & ,rw(1,1,nz+1-kts(3)),rw(1,1,nz+1-kts(4)),rw(1,1,nz+1-kts(5))
c     & ,rw(1,1,nz+1-kts(6)),rw(1,1,nz+1-kts(7)),rw(1,1,nz+1-kts(8))
c     & ,rw(1,1,nz+1-kts(9)),rw(1,1,nz+1-kts(10)),rw(1,1,nz+1-kts(11))
c      write(80,300)rtime,rv(nx/2,ny/2,kts(1)),rv(nx/2,ny/2,kts(2))
c     & ,rv(nx/2,ny/2,kts(3)),rv(nx/2,ny/2,kts(4)),rv(nx/2,ny/2,kts(5))
c     & ,rv(nx/2,ny/2,kts(6)),rv(nx/2,ny/2,kts(7)),rv(nx/2,ny/2,kts(8))
c     & ,rv(nx/2,ny/2,kts(9)),rv(nx/2,ny/2,kts(10)),rv(nx/2,ny/2,kts(11))
c      write(81,300)rtime,ru(nx/2,ny/2,kts(1)),ru(nx/2,ny/2,kts(2))
c     & ,ru(nx/2,ny/2,kts(3)),ru(nx/2,ny/2,kts(4)),ru(nx/2,ny/2,kts(5))
c     & ,ru(nx/2,ny/2,kts(6)),ru(nx/2,ny/2,kts(7)),ru(nx/2,ny/2,kts(8))
c     & ,ru(nx/2,ny/2,kts(9)),ru(nx/2,ny/2,kts(10)),ru(nx/2,ny/2,kts(11))
cC$OMP SECTION
c      write(82,300)rtime,rw(nx/2,ny/2,kts(1)),rw(nx/2,ny/2,kts(2))
c     & ,rw(nx/2,ny/2,kts(3)),rw(nx/2,ny/2,kts(4)),rw(nx/2,ny/2,kts(5))
c     & ,rw(nx/2,ny/2,kts(6)),rw(nx/2,ny/2,kts(7)),rw(nx/2,ny/2,kts(8))
c     & ,rw(nx/2,ny/2,kts(9)),rw(nx/2,ny/2,kts(10)),rw(nx/2,ny/2,kts(11))
c      write(90,300)rtime,rv(1,ny/2,nz+1-kts(1)),rv(1,ny/2,nz+1-kts(2))
c     & ,rv(1,ny/2,nz+1-kts(3)),rv(1,ny/2,nz+1-kts(4))
c     & ,rv(1,ny/2,nz+1-kts(5)),rv(1,ny/2,nz+1-kts(6))
c     & ,rv(1,ny/2,nz+1-kts(7)),rv(1,ny/2,nz+1-kts(8))
c     & ,rv(1,ny/2,nz+1-kts(9)),rv(1,ny/2,nz+1-kts(10))
c     & ,rv(1,ny/2,nz+1-kts(11))
c      write(91,300)rtime,ru(1,ny/2,nz+1-kts(1)),ru(1,ny/2,nz+1-kts(2))
c     & ,ru(1,ny/2,nz+1-kts(3)),ru(1,ny/2,nz+1-kts(4))
c     & ,ru(1,ny/2,nz+1-kts(5)),ru(1,ny/2,nz+1-kts(6))
c     & ,ru(1,ny/2,nz+1-kts(7)),ru(1,ny/2,nz+1-kts(8))
c     & ,ru(1,ny/2,nz+1-kts(9)),ru(1,ny/2,nz+1-kts(10))
c     & ,ru(1,ny/2,nz+1-kts(11))
c      write(92,300)rtime,rw(1,ny/2,nz+1-kts(1)),rw(1,ny/2,nz+1-kts(2))
c     & ,rw(1,ny/2,nz+1-kts(3)),rw(1,ny/2,nz+1-kts(4))
c     & ,rw(1,ny/2,nz+1-kts(5)),rw(1,ny/2,nz+1-kts(6))
c     & ,rw(1,ny/2,nz+1-kts(7)),rw(1,ny/2,nz+1-kts(8))
c     & ,rw(1,ny/2,nz+1-kts(9)),rw(1,ny/2,nz+1-kts(10))
c     & ,rw(1,ny/2,nz+1-kts(11))
cC$OMP SECTION
c      if(FLTHRM)then
c       write(63,300)rtime,rtt(nx/2,1,kts(1)),rtt(nx/2,1,kts(2))
c     &  ,rtt(nx/2,1,kts(3)),rtt(nx/2,1,kts(4)),rtt(nx/2,1,kts(5))
c     &  ,rtt(nx/2,1,kts(6)),rtt(nx/2,1,kts(7)),rtt(nx/2,1,kts(8))
c     &  ,rtt(nx/2,1,kts(9)),rtt(nx/2,1,kts(10)),rtt(nx/2,1,kts(11))
c       write(73,300)rtime,rtt(1,1,nz+1-kts(1)),rtt(1,1,nz+1-kts(2))
c     &  ,rtt(1,1,nz+1-kts(3)),rtt(1,1,nz+1-kts(4)),rtt(1,1,nz+1-kts(5))
c     &  ,rtt(1,1,nz+1-kts(6)),rtt(1,1,nz+1-kts(7)),rtt(1,1,nz+1-kts(8))
c     &  ,rtt(1,1,nz+1-kts(9)),rtt(1,1,nz+1-kts(10))
c     &  ,rtt(1,1,nz+1-kts(11))
c       write(83,300)rtime,rtt(nx/2,ny/2,kts(1)),rtt(nx/2,ny/2,kts(2))
c     &  ,rtt(nx/2,ny/2,kts(3)),rtt(nx/2,ny/2,kts(4))
c     &  ,rtt(nx/2,ny/2,kts(5)),rtt(nx/2,ny/2,kts(6))
c     &  ,rtt(nx/2,ny/2,kts(7)),rtt(nx/2,ny/2,kts(8))
c     &  ,rtt(nx/2,ny/2,kts(9)),rtt(nx/2,ny/2,kts(10))
c     &  ,rtt(nx/2,ny/2,kts(11))
c       write(93,300)rtime
c     &  ,rtt(1,ny/2,nz+1-kts(1)),rtt(1,ny/2,nz+1-kts(2))
c     &  ,rtt(1,ny/2,nz+1-kts(3)),rtt(1,ny/2,nz+1-kts(4))
c     &  ,rtt(1,ny/2,nz+1-kts(5)),rtt(1,ny/2,nz+1-kts(6))
c     &  ,rtt(1,ny/2,nz+1-kts(7)),rtt(1,ny/2,nz+1-kts(8))
c     &  ,rtt(1,ny/2,nz+1-kts(9)),rtt(1,ny/2,nz+1-kts(10))
c     &  ,rtt(1,ny/2,nz+1-kts(11))
c      endif
cC$OMP END SECTIONS NOWAIT
cC$OMP END PARALLEL
c
c 300  format(g25.15,g25.15,g25.15,g25.15,g25.15,g25.15,g25.15,g25.15
c     &,g25.15,g25.15,g25.15,g25.15)
 300  format(12(3x,g15.8))
c
c     decide weather accumulate 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc     1st version
cc     this version does not work well for restart and for np=1
c      if(mod(itime,dnp).ne.0)return
cc
cc     update counter
c      if(mod(itime/dnp,np).eq.1)NT=NT+1
cc
cc     set the phase index where to accumulate 
c      index=itime
c      if(itime.gt.np*dnp)index=itime-(NT-1)*np*dnp
c      ip=index/dnp
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc     2nd version
cc     this version does not work well for restart
c      if(np.eq.1)then
c         NT=NT+1
c         ip=1
c      else
c         if(mod(itime,dnp).ne.0)return
cc     update counter
c         if(mod(itime/dnp,np).eq.1)NT=NT+1
cc     set the phase index where to accumulate 
c         index=itime
c         if(itime.gt.np*dnp)index=itime-(NT-1)*np*dnp
c         ip=index/dnp
c         if(ip.gt.np.or.ip.lt.1)then
c            write(26,*)'Wrong phase in make_means'
c            write(*,*)'Wrong phase in make_means'
c         endif
c      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     3rd version slves np=1 and restart at any interation
      if(mod(itime,dnp).ne.0)return
      itimeaux=mod(itime,np*dnp)
      itimeaux=itimeaux/dnp+1
      ip=mod(itimeaux+(np-2),np)+1
c
c no need to set special case of np=1, previous algorithm works
c      if(np.eq.1)ip=1
c     
      if(ip.gt.np.or.ip.lt.1)then
         write(26,*)'Wrong phase in make_means'
         write(26,*)'np:',np,'  ip:',ip,'  NT:',NT(ip)
         write(26,*)'dnp:',dnp,'  itime:',itime
         write( *,*)'Wrong phase in make_means'
         write( *,*)'np:',np,'  ip:',ip,'  NT:',NT(ip)
         write( *,*)'dnp:',dnp,'  itime:',itime
         stop
      endif
c
      NT(ip)=NT(ip)+1
c
c      write(17,*)'itime:',itime,'  ip:',ip,'  NT:',NT(ip)
c      write( *,*)'itime:',itime,'  ip:',ip,'  NT:',NT(ip)
c  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     output adcp data here for Jr
c      call output_adcpdata
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     means
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sum,svm,swm,su2m,sv2m,sw2m,
C$OMP&               su3m,sv3m,sw3m,su4m,sv4m,sw4m,
C$OMP&               suvm,suwm,svwm,
C$OMP&               spm,spwm,su2wm,sv2wm)
      do k=1,nz
         sum   = 0.d0 
         svm   = 0.d0
         swm   = 0.d0
         su2m  = 0.d0
         sv2m  = 0.d0
         sw2m  = 0.d0
         su3m  = 0.d0
         sv3m  = 0.d0
         sw3m  = 0.d0
         su4m  = 0.d0
         sv4m  = 0.d0
         sw4m  = 0.d0
         suvm  = 0.d0
         suwm  = 0.d0
         svwm  = 0.d0
         spm   = 0.d0
         spwm  = 0.d0
         su2wm = 0.d0
         sv2wm = 0.d0
         do j=1,ny
            do i=1,nx
               sum   = sum   + ru( i,j,k) 
               svm   = svm   + rv( i,j,k)
               swm   = swm   + rw( i,j,k)
               su2m  = su2m  + ru( i,j,k)**2
               sv2m  = sv2m  + rv( i,j,k)**2
               sw2m  = sw2m  + rw( i,j,k)**2
               su3m  = su3m  + ru( i,j,k)**3
               sv3m  = sv3m  + rv( i,j,k)**3
               sw3m  = sw3m  + rw( i,j,k)**3
               su4m  = su4m  + ru( i,j,k)**4
               sv4m  = sv4m  + rv( i,j,k)**4
               sw4m  = sw4m  + rw( i,j,k)**4
               suvm  = suvm  + ru( i,j,k)*rv( i,j,k)
               suwm  = suwm  + ru( i,j,k)*rw( i,j,k)
               svwm  = svwm  + rv( i,j,k)*rw( i,j,k)
               spm   = spm   + rp_dum( i,j,k)
               spwm  = spwm  + rp_dum( i,j,k)*rw( i,j,k)
               su2wm = su2wm + ru( i,j,k)**2*rw( i,j,k)
               sv2wm = sv2wm + rv( i,j,k)**2*rw( i,j,k)
            enddo
         enddo
         um(  k,ip) = um(  k,ip) + sum   / dble(nxny) 
         vm(  k,ip) = vm(  k,ip) + svm   / dble(nxny) 
         wm(  k,ip) = wm(  k,ip) + swm   / dble(nxny) 
         u2m( k,ip) = u2m( k,ip) + su2m  / dble(nxny)
         v2m( k,ip) = v2m( k,ip) + sv2m  / dble(nxny)
         w2m( k,ip) = w2m( k,ip) + sw2m  / dble(nxny)
         u3m( k,ip) = u3m( k,ip) + su3m  / dble(nxny)
         v3m( k,ip) = v3m( k,ip) + sv3m  / dble(nxny)
         w3m( k,ip) = w3m( k,ip) + sw3m  / dble(nxny)
         u4m( k,ip) = u4m( k,ip) + su4m  / dble(nxny)
         v4m( k,ip) = v4m( k,ip) + sv4m  / dble(nxny)
         w4m( k,ip) = w4m( k,ip) + sw4m  / dble(nxny)
         uvm( k,ip) = uvm( k,ip) + suvm  / dble(nxny)
         uwm( k,ip) = uwm( k,ip) + suwm  / dble(nxny)
         vwm( k,ip) = vwm( k,ip) + svwm  / dble(nxny)
         pm(  k,ip) = pm(  k,ip) + spm   / dble(nxny)
         pwm( k,ip) = pwm( k,ip) + spwm  / dble(nxny)
         u2wm(k,ip) = u2wm(k,ip) + su2wm / dble(nxny)
         v2wm(k,ip) = v2wm(k,ip) + sv2wm / dble(nxny)
      enddo
C$OMP END PARALLEL DO
c
      if(FLTHRM)then
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sttm,stt2m,stt3m,stt4m,suttm,svttm,swttm)
         do k=1,nz
            sttm  = 0.d0
            stt2m = 0.d0
            stt3m = 0.d0
            stt4m = 0.d0
            suttm = 0.d0
            svttm = 0.d0
            swttm = 0.d0 
            do j=1,ny
               do i=1,nx
                  sttm  = sttm  + rtt(i,j,k)
                  stt2m = stt2m + rtt(i,j,k)**2
                  stt3m = stt3m + rtt(i,j,k)**3
                  stt4m = stt4m + rtt(i,j,k)**4
                  suttm = suttm + ru( i,j,k)*rtt(i,j,k)
                  svttm = svttm + rv( i,j,k)*rtt(i,j,k)
                  swttm = swttm + rw( i,j,k)*rtt(i,j,k)
               enddo
            enddo
            ttm( k,ip) = ttm( k,ip) + sttm  / dble(nxny) 
            tt2m(k,ip) = tt2m(k,ip) + stt2m / dble(nxny)
            tt3m(k,ip) = tt3m(k,ip) + stt3m / dble(nxny)
            tt4m(k,ip) = tt4m(k,ip) + stt4m / dble(nxny)
            uttm(k,ip) = uttm(k,ip) + suttm / dble(nxny)
            vttm(k,ip) = vttm(k,ip) + svttm / dble(nxny)
            wttm(k,ip) = wttm(k,ip) + swttm / dble(nxny)
         enddo
C$OMP END PARALLEL DO
      endif
c
c     output some check-up data
      aux =um( 3*nz/4,ip)/dble(NT(ip))
      avx =vm( 3*nz/4,ip)/dble(NT(ip))
      awx =wm( 3*nz/4,ip)/dble(NT(ip))
      attx=ttm(3*nz/4,ip)/dble(NT(ip))
      write(10,10) rtime,aux,avx,awx,attx
     &            ,u2m( 3*nz/4,ip)/dble(NT(ip))- aux**2
     &            ,v2m( 3*nz/4,ip)/dble(NT(ip))- avx**2
     &            ,w2m( 3*nz/4,ip)/dble(NT(ip))- awx**2
     &            ,tt2m(3*nz/4,ip)/dble(NT(ip))-attx**2
c 10   format(g25.15,g25.15,g25.15,g25.15,g25.15,g25.15,g25.15,g25.15
c     &      ,g25.15)
 10   format(9(3x,g15.8))
c     
c     z-component of gradients
      call ddz(ru,rp_dum)
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,suuzm)
      do k=1,nz
         suuzm = 0.d0 
         do j=1,ny
            do i=1,nx
               suuzm = suuzm + ru(i,j,k)*rp_dum(i,j,k) 
            enddo
         enddo
         uuzm(k,ip)  = uuzm(k,ip) + suuzm / dble(nxny) 
      enddo
C$OMP END PARALLEL DO
      call ddz(rv,ru)
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,svvzm)
      do k=1,nz
         svvzm = 0.d0 
         do j=1,ny
            do i=1,nx
               svvzm = svvzm + rv(i,j,k)*ru(i,j,k) 
            enddo
         enddo
         vvzm(k,ip)  = vvzm(k,ip) + svvzm / dble(nxny) 
      enddo
C$OMP END PARALLEL DO
      call ddz(rw,rv)
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,swwzm)
      do k=1,nz
         swwzm = 0.d0 
         do j=1,ny
            do i=1,nx
               swwzm = swwzm + rw(i,j,k)*rv(i,j,k) 
            enddo
         enddo
         wwzm(k,ip)  = wwzm(k,ip) + swwzm / dble(nxny) 
      enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sum,svm,swm,su2m,sv2m,sw2m,
C$OMP&               su3m,sv3m,sw3m,su4m,sv4m,sw4m)
      do k=1,nz
         sum   = 0.d0 
         su2m  = 0.d0
         su3m  = 0.d0
         su4m  = 0.d0
         svm   = 0.d0
         sv2m  = 0.d0
         sv3m  = 0.d0
         sv4m  = 0.d0
         swm   = 0.d0
         sw2m  = 0.d0
         sw3m  = 0.d0
         sw4m  = 0.d0
         do j=1,ny
            do i=1,nx
               sum   = sum   + rp_dum( i,j,k) 
               su2m  = su2m  + rp_dum( i,j,k)**2
               su3m  = su3m  + rp_dum( i,j,k)**3
               su4m  = su4m  + rp_dum( i,j,k)**4
c
               svm   = svm   + ru( i,j,k)
               sv2m  = sv2m  + ru( i,j,k)**2
               sv3m  = sv3m  + ru( i,j,k)**3
               sv4m  = sv4m  + ru( i,j,k)**4
c     
               swm   = swm   + rv( i,j,k)
               sw2m  = sw2m  + rv( i,j,k)**2
               sw3m  = sw3m  + rv( i,j,k)**3
               sw4m  = sw4m  + rv( i,j,k)**4
            enddo
         enddo
         uzm(  k,ip) = uzm(  k,ip) + sum   / dble(nxny) 
         uz2m( k,ip) = uz2m( k,ip) + su2m  / dble(nxny)
         uz3m( k,ip) = uz3m( k,ip) + su3m  / dble(nxny)
         uz4m( k,ip) = uz4m( k,ip) + su4m  / dble(nxny)
c     
         vzm(  k,ip) = vzm(  k,ip) + svm   / dble(nxny) 
         vz2m( k,ip) = vz2m( k,ip) + sv2m  / dble(nxny)
         vz3m( k,ip) = vz3m( k,ip) + sv3m  / dble(nxny)
         vz4m( k,ip) = vz4m( k,ip) + sv4m  / dble(nxny)
c     
         wzm(  k,ip) = wzm(  k,ip) + swm   / dble(nxny) 
         wz2m( k,ip) = wz2m( k,ip) + sw2m  / dble(nxny)
         wz3m( k,ip) = wz3m( k,ip) + sw3m  / dble(nxny)
         wz4m( k,ip) = wz4m( k,ip) + sw4m  / dble(nxny)
      enddo
C$OMP END PARALLEL DO
c     
      if(FLTHRM)then
         call ddz(rtt,rw)
c
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sttm,stt2m,stt3m,stt4m)
         do k=1,nz
            sttm  = 0.d0
            stt2m = 0.d0
            stt3m = 0.d0
            stt4m = 0.d0
            do j=1,ny
               do i=1,nx
                  sttm  = sttm  + rw(i,j,k)
                  stt2m = stt2m + rw(i,j,k)**2
                  stt3m = stt3m + rw(i,j,k)**3
                  stt4m = stt4m + rw(i,j,k)**4
               enddo
            enddo
            ttzm( k,ip) = ttzm( k,ip) + sttm  / dble(nxny) 
            ttz2m(k,ip) = ttz2m(k,ip) + stt2m / dble(nxny)
            ttz3m(k,ip) = ttz3m(k,ip) + stt3m / dble(nxny)
            ttz4m(k,ip) = ttz4m(k,ip) + stt4m / dble(nxny)            
         enddo
C$OMP END PARALLEL DO
      endif
c
c     output u*^2 as check-up data
      write(11,11) rtime,vzm(nz,ip)/dble(NT(ip))/re
     &                  ,vzm( 1,ip)/dble(NT(ip))/re
     &                  ,uzm(nz,ip)/dble(NT(ip))/re
     &                  ,uzm( 1,ip)/dble(NT(ip))/re
c 11   format(g25.15,g25.15,g25.15,g25.15,g25.15)
 11   format(5(3x,g15.8))
c
c     x-component of gradients
      call ddx_f_new(u,ru)
      call ddx_f_new(v,rv)
      call ddx_f_new(w,rw)
c     
      call fft_f2r_2d_new(ru)
      call fft_f2r_2d_new(rv)
      call fft_f2r_2d_new(rw)
c
c     produce u again needed here
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k)
      do k=1,nz
         call dcopy(nx0y0,u(1,1,k),1,rp_dum(1,1,k),1)
      enddo
C$OMP END PARALLEL DO
      call fft_f2r_2d_new(rp_dum)
c
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sum,svm,swm,su2m,sv2m,sw2m,
C$OMP&               su3m,sv3m,sw3m,su4m,sv4m,sw4m,
C$OMP&               suwxm)
      do k=1,nz
         sum   = 0.d0 
         su2m  = 0.d0
         su3m  = 0.d0
         su4m  = 0.d0
         svm   = 0.d0
         sv2m  = 0.d0
         sv3m  = 0.d0
         sv4m  = 0.d0
         swm   = 0.d0
         sw2m  = 0.d0
         sw3m  = 0.d0
         sw4m  = 0.d0
         suwxm = 0.d0 
         do j=1,ny
            do i=1,nx
               sum   = sum   + ru( i,j,k) 
               su2m  = su2m  + ru( i,j,k)**2
               su3m  = su3m  + ru( i,j,k)**3
               su4m  = su4m  + ru( i,j,k)**4
c     
               svm   = svm   + rv( i,j,k)
               sv2m  = sv2m  + rv( i,j,k)**2
               sv3m  = sv3m  + rv( i,j,k)**3
               sv4m  = sv4m  + rv( i,j,k)**4
c     
               swm   = swm   + rw( i,j,k)
               sw2m  = sw2m  + rw( i,j,k)**2
               sw3m  = sw3m  + rw( i,j,k)**3
               sw4m  = sw4m  + rw( i,j,k)**4
c
               suwxm = suwxm + rp_dum( i,j,k)*rw( i,j,k) 
            enddo
         enddo
         uxm(  k,ip) = uxm(  k,ip) + sum   / dble(nxny) 
         ux2m( k,ip) = ux2m( k,ip) + su2m  / dble(nxny)
         ux3m( k,ip) = ux3m( k,ip) + su3m  / dble(nxny)
         ux4m( k,ip) = ux4m( k,ip) + su4m  / dble(nxny)
c     
         vxm(  k,ip) = vxm(  k,ip) + svm   / dble(nxny) 
         vx2m( k,ip) = vx2m( k,ip) + sv2m  / dble(nxny)
         vx3m( k,ip) = vx3m( k,ip) + sv3m  / dble(nxny)
         vx4m( k,ip) = vx4m( k,ip) + sv4m  / dble(nxny)
c     
         wxm(  k,ip) = wxm(  k,ip) + swm   / dble(nxny) 
         wx2m( k,ip) = wx2m( k,ip) + sw2m  / dble(nxny)
         wx3m( k,ip) = wx3m( k,ip) + sw3m  / dble(nxny)
         wx4m( k,ip) = wx4m( k,ip) + sw4m  / dble(nxny)
c
         uwxm(k,ip)  = uwxm( k,ip) + suwxm / dble(nxny) 
      enddo
C$OMP END PARALLEL DO
c
      if(FLTHRM)then 
         call ddx_f_new(tt,rtt)
         call fft_f2r_2d_new(rtt)
c
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sttm,stt2m,stt3m,stt4m)
         do k=1,nz
            sttm  = 0.d0
            stt2m = 0.d0
            stt3m = 0.d0
            stt4m = 0.d0
            do j=1,ny
               do i=1,nx
                  sttm  = sttm  + rtt(i,j,k)
                  stt2m = stt2m + rtt(i,j,k)**2
                  stt3m = stt3m + rtt(i,j,k)**3
                  stt4m = stt4m + rtt(i,j,k)**4
               enddo
            enddo
            ttxm( k,ip) = ttxm( k,ip) + sttm  / dble(nxny) 
            ttx2m(k,ip) = ttx2m(k,ip) + stt2m / dble(nxny)
            ttx3m(k,ip) = ttx3m(k,ip) + stt3m / dble(nxny)
            ttx4m(k,ip) = ttx4m(k,ip) + stt4m / dble(nxny)
         enddo
C$OMP END PARALLEL DO
      endif
c
c     y-component of gradients
      call ddy_f_new(u,ru)
      call ddy_f_new(v,rv)
      call ddy_f_new(w,rw)
c     
      call fft_f2r_2d_new(ru)
      call fft_f2r_2d_new(rv)
      call fft_f2r_2d_new(rw)
c
c     produce v again needed here
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k)
      do k=1,nz
         call dcopy(nx0y0,v(1,1,k),1,rp_dum(1,1,k),1)
      enddo
C$OMP END PARALLEL DO
      call fft_f2r_2d_new(rp_dum)
c
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sum,svm,swm,su2m,sv2m,sw2m,
C$OMP&               su3m,sv3m,sw3m,su4m,sv4m,sw4m,
C$OMP&               svwym)
      do k=1,nz
         sum   = 0.d0 
         su2m  = 0.d0
         su3m  = 0.d0
         su4m  = 0.d0
         svm   = 0.d0
         sv2m  = 0.d0
         sv3m  = 0.d0
         sv4m  = 0.d0
         swm   = 0.d0
         sw2m  = 0.d0
         sw3m  = 0.d0
         sw4m  = 0.d0
         svwym = 0.d0
         do j=1,ny
            do i=1,nx
               sum   = sum   + ru( i,j,k) 
               su2m  = su2m  + ru( i,j,k)**2
               su3m  = su3m  + ru( i,j,k)**3
               su4m  = su4m  + ru( i,j,k)**4
c     
               svm   = svm   + rv( i,j,k)
               sv2m  = sv2m  + rv( i,j,k)**2
               sv3m  = sv3m  + rv( i,j,k)**3
               sv4m  = sv4m  + rv( i,j,k)**4
c     
               swm   = swm   + rw( i,j,k)
               sw2m  = sw2m  + rw( i,j,k)**2
               sw3m  = sw3m  + rw( i,j,k)**3
               sw4m  = sw4m  + rw( i,j,k)**4
c     
               svwym = svwym + rp_dum( i,j,k)*rw( i,j,k)
            enddo
         enddo
         uym(  k,ip) = uym(  k,ip) + sum   / dble(nxny) 
         uy2m( k,ip) = uy2m( k,ip) + su2m  / dble(nxny)
         uy3m( k,ip) = uy3m( k,ip) + su3m  / dble(nxny)
         uy4m( k,ip) = uy4m( k,ip) + su4m  / dble(nxny)
c     
         vym(  k,ip) = vym(  k,ip) + svm   / dble(nxny) 
         vy2m( k,ip) = vy2m( k,ip) + sv2m  / dble(nxny)
         vy3m( k,ip) = vy3m( k,ip) + sv3m  / dble(nxny)
         vy4m( k,ip) = vy4m( k,ip) + sv4m  / dble(nxny)
c     
         wym(  k,ip) = wym(  k,ip) + swm   / dble(nxny) 
         wy2m( k,ip) = wy2m( k,ip) + sw2m  / dble(nxny)
         wy3m( k,ip) = wy3m( k,ip) + sw3m  / dble(nxny)
         wy4m( k,ip) = wy4m( k,ip) + sw4m  / dble(nxny)
c     
         vwym( k,ip) = vwym( k,ip) + svwym / dble(nxny) 
      enddo
C$OMP END PARALLEL DO
c
      if(FLTHRM)then 
         call ddy_f_new(tt,rtt)
         call fft_f2r_2d_new(rtt)
c
C$OMP PARALLEL DO DEFAULT(SHARED),
C$OMP& PRIVATE(i,j,k,sttm,stt2m,stt3m,stt4m)
         do k=1,nz
            sttm  = 0.d0
            stt2m = 0.d0
            stt3m = 0.d0
            stt4m = 0.d0
            do j=1,ny
               do i=1,nx
                  sttm  = sttm  + rtt(i,j,k)
                  stt2m = stt2m + rtt(i,j,k)**2
                  stt3m = stt3m + rtt(i,j,k)**3
                  stt4m = stt4m + rtt(i,j,k)**4
               enddo
            enddo
            ttym( k,ip) = ttym( k,ip) + sttm  / dble(nxny) 
            tty2m(k,ip) = tty2m(k,ip) + stt2m / dble(nxny)
            tty3m(k,ip) = tty3m(k,ip) + stt3m / dble(nxny)
            tty4m(k,ip) = tty4m(k,ip) + stt4m / dble(nxny)
         enddo
C$OMP END PARALLEL DO
      endif
c
c     u-spectra
      call make_spec_new(u,espx,dspx,espy,dspy,xrms,txrms,
     &     yrms,tyrms)
c     
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
         do i=1,nxh
            uespxm(i,k,ip)= uespxm(i,k,ip)+ espx(i,k)
c            udspxm(i,k,ip)= udspxm(i,k,ip)+ dspx(i,k)
         enddo
         do j=1,nyh
            uespym(j,k,ip)= uespym(j,k,ip)+ espy(j,k)
c            udspym(j,k,ip)= udspym(j,k,ip)+ dspy(j,k)
         enddo
      enddo
C$OMP END PARALLEL DO
c     
c     v-spectra
      call make_spec_new(v,espx,dspx,espy,dspy,xrms,txrms,
     &     yrms,tyrms)
c     
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
         do i=1,nxh
            vespxm(i,k,ip)= vespxm(i,k,ip)+ espx(i,k)
c            vdspxm(i,k,ip)= vdspxm(i,k,ip)+ dspx(i,k)
         enddo
         do j=1,nyh
            vespym(j,k,ip)= vespym(j,k,ip)+ espy(j,k)
c            vdspym(j,k,ip)= vdspym(j,k,ip)+ dspy(j,k)
         enddo
      enddo
C$OMP END PARALLEL DO
c
c     w-spectra
      call make_spec_new(w,espx,dspx,espy,dspy,xrms,txrms,
     &     yrms,tyrms)
c     
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k)
      do k=1,nz
         do i=1,nxh
            wespxm(i,k,ip)= wespxm(i,k,ip)+ espx(i,k)
c            wdspxm(i,k,ip)= wdspxm(i,k,ip)+ dspx(i,k)
         enddo
         do j=1,nyh
            wespym(j,k,ip)= wespym(j,k,ip)+ espy(j,k)
c            wdspym(j,k,ip)= wdspym(j,k,ip)+ dspy(j,k)
         enddo
      enddo
C$OMP END PARALLEL DO
c
      return
      END
