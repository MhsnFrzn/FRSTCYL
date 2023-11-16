      PROGRAM main
#include "common.inc"
c
      INTEGER ch_fin,m_fin,sph_fin,itfirst,iend,i,j,k
      INTEGER tt_yes,tt2_yes,ibm_yes,mov_yes,lag_yes,ngbal,n1d,code
      REAL*8  dt0,dt_p,KE,Ep,diss,ddtEp
      REAL*8  espz(nz) 
      CHARACTER*32 identifier,name_dudn,name_height,itime_str 
c
      REAL*8 espx(nxh,nz),dspx(nxh,nz),espy(nyh,nz),dspy(nyh,nz)
      REAL*8 fmean(nz),xrms(nz),yrms(nz),pmean(nz),prms(nz)
      REAL*8 tfmean,txrms,tyrms,tpmean,tprms !,igmod
      REAL*8 auxdepo(nx0,ny0),auxeros(nx0,ny0)
      REAL*8 u_dum(nx0,ny0,nz),v_dum(nx0,ny0,nz),w_dum(nx0,ny0,nz)

      print *, 'Program running'
c
c     read in the parameters for the time-stepping and the filenames
c     ch_fin = 0 generates an initial solution
c     ch_fin = 1 reads the initial data from ch_file.(istart)
c
c     outfirst = whether to output data (other than c and u) at
c     first time-step
c
      read '(a)', identifier
      if (identifier(1:5).ne.'#MAIN') then
        print *, 'Expecting identifier #MAIN in input file'
        stop
      endif
      read '(a)', ch_file 
      read(*,*) rlenx,rleny,rlenz,uvbb,uvbt
      read(*,*) istart,iters,ch_fin,m_fin
      read(*,*) oi_chan,out_press
      read(*,*) oi_timer,oi_cfl,oi_spec,oi_mean,oi_gbal,oi_1d
      read(*,*) re,ppA,dt_p,cflmax,ppT,theta
      read(*,*) kts(1),kts(2),kts(3),kts(4),kts(5),kts(6),kts(7),kts(8)
     &         ,kts(9),kts(10),kts(11)
      read(*,*) ibm_yes
      read(*,*) mov_yes
      read(*,*) lag_yes
      read(*,*) tt_yes
      read(*,*) tt2_yes
      read(*,*) AF
      read(*,*) diff_coef
c
      dt = dt_p
      FLTIMER = oi_timer.gt.0
      FLAVER  = oi_mean.gt.0
      FLIBM   = ibm_yes.ne.0
      FLMOV   = mov_yes.ne.0
      FLSURF  = lag_yes.ne.0
      FLTHRM  = tt_yes.ne.0
      FLTHRM2 = tt2_yes.ne.0
      FLPGVAR = ppT.gt.0.d0
c
      if (FLTIMER) then
        call init_timers
      else
        oi_timer = 1
      endif
c
      read '(a)', identifier
      if (identifier(1:9).ne.'#END_MAIN') then
        print *, 'Expecting identifier #END_MAIN in input file'
        stop
      endif
c
      if (FLIBM) then
         read '(a)', identifier
         if (identifier(1:5).ne.'#SURF') then
            print *, 'Expecting identifier #SURF in input file'
            stop
         endif
         read '(a)', lag_file
         read(*,*) T_wall, T_inf, Y_wall, Y_inf!, eps_bed
         read(*,*) zplus, rplus
         read '(a)', identifier
         if (identifier(1:9).ne.'#END_SURF') then
            print *, 'Expecting identifier #END_SURF in input file'
            stop
         endif
      endif

c
      if (FLTHRM) then
         read '(a)', identifier
         if (identifier(1:5).ne.'#THRM') then
            print *, 'Expecting identifier #THRM in input file'
            stop
         endif
         read '(a)', tt_file
         read(*,*) pr,ra,Tmax
         read(*,*) ttab(1),ttat(1) 
         read(*,*) ttbb(1),ttbt(1)
         read(*,*) ttgb(1),ttgt(1)
         read(*,*) ws,rey_p
         read '(a)', identifier
         if (identifier(1:9).ne.'#END_THRM') then
            print *, 'Expecting identifier #END_THRM in input file'
            stop
         endif
      endif
c
      if (FLTHRM2) then
         read '(a)', identifier
         if (identifier(1:6).ne.'#THRM2') then
            print *, 'Expecting identifier #THRM2 in input file'
            stop
         endif
         read '(a)', tt2_file
         read(*,*) pr2,ra2,Tmax2
         read(*,*) tt2ab(1),tt2at(1)
         read(*,*) tt2bb(1),tt2bt(1)
         read(*,*) tt2gb(1),tt2gt(1)
c         read(*,*) ws,rey_p
         read '(a)', identifier
         if (identifier(1:10).ne.'#END_THRM2') then
            print *, 'Expecting identifier #END_THRM2 in input file'
            stop
         endif
      endif

c
      FLSTLE  = ws.ne.0.0
c
      epssvx=0.d0
      epssvy=0.d0
      if (FLTHRM) then
         read '(a)', identifier
         if (identifier(1:8).ne.'#SP_VISC') then
            print *, 'Expecting identifier #SP_VISC in input file'
            stop
         endif
         read(*,*) epssvx,epssvy
         read(*,*) msvx,msvy
         read '(a)', identifier
         if (identifier(1:12).ne.'#END_SP_VISC') then
            print *, 'Expecting identifier #END_SP_VISC in input file'
            stop
         endif
      endif
      FLSV_NO=.false.
      FLSV_YES=epssvx.gt.0.or.epssvy.gt.0
c
      epssvx2=0.d0
      epssvy2=0.d0     
      if (FLTHRM2) then
         read '(a)', identifier
         if (identifier(1:9).ne.'#SP_VISC2') then
            print *, 'Expecting identifier #SP_VISC2 in input file'
            stop
         endif
         read(*,*) epssvx2,epssvy2
         read(*,*) msvx2,msvy2 
         read '(a)', identifier
         if (identifier(1:13).ne.'#END_SP_VISC2') then
            print *, 'Expecting identifier #END_SP_VISC2 in input file'
            stop
         endif
      endif         
      FLSV2_NO=.false.
      FLSV2_YES=epssvx2.gt.0.or.epssvy2.gt.0
c
      print *, 'Input read'
c
      itfirst = istart+1
      iend = istart + iters
      itime = istart
c
      print *,              'Flow file root:           ', ch_file
      if (FLTHRM ) print *, 'Temperature file root:    ', tt_file
      if (FLTHRM2) print *, 'Humidity    file root:    ', tt2_file
c
      open(26,file='logfile')
c      open(31,file='gbal')
c      rewind(31)
c      ngbal = 0
c      open(30,file='1dout')
c      rewind(30)
c      n1d = 0
c
      
      call initial(dt0,ch_fin,m_fin)
c
c      call comp_flux
      
c      stop
c produce complete output for initial condition if ...
      if(ch_fin.eq.0)then
c         call make_gbalrhs(KE,Ep,diss,ddtEp)
c         write(31,200) rtime,KE,Ep,diss,ddtEp
c         ngbal = ngbal + 1
c     
         call output
c     
         if(FLAVER)then
         endif
      endif !(ch_fin.eq.0)
c     
      call divg
c
c Nadim need to include info on tt2 in the logfile below
      write(26,*)' rlenx, rleny, rlenz, dt'
      write(26,*)  rlenx, rleny, rlenz, dt
      write(26,*)' deltax, deltay, dt0'
      write(26,*)  deltax, deltay, dt0
      write(26,*)' uvab, uvat, uvbb, uvbt'
      write(26,*)  uvab, uvat, uvbb, uvbt
      write(26,*)' istart, iters, itfirst, iend'
      write(26,*)  istart, iters, itfirst, iend
      write(26,*)' oi_chan, out_press'
      write(26,*)  oi_chan, out_press
      write(26,*)' oi_timer, oi_cfl, oi_spec, oi_mean, oi_gbal, oi_1d'
      write(26,*)  oi_timer, oi_cfl, oi_spec, oi_mean, oi_gbal, oi_1d
      write(26,*)' nx, ny, nz, nxl, nyl'
      write(26,*)  nx, ny, nz, nxl, nyl
c      write(26,*)' Re,gradpbar,dt_p,cflmax'
      write(26,*)' Re,ppA,dt_p,cflmax,ppT,theta'
c      write(26,*)  re,gradpbar,dt_p,cflmax
      write(26,*)  re,ppA,dt_p,cflmax,ppT,theta
      write(26,*) 'kts(1) to kts(11)'
      write(26,*) kts(1),kts(2),kts(3),kts(4),kts(5),kts(6)
     &           ,kts(7),kts(8),kts(9),kts(10),kts(11)
      write(26,*)' tt_yes'
      write(26,*)  tt_yes
      if (FLTHRM) then
        write(26,*)' Pr, Ra, Tmax'
        write(26,*)  pr, ra, Tmax
        write(26,*)' BC: ttat(1), ttab(1), ttbt(1), ttbb(1)
     &, ttgt(1), ttgb(1)'
        write(26,*)      ttat(1), ttab(1), ttbt(1), ttbb(1)
     &, ttgt(1), ttgb(1)
        write(26,*)' ws, rey_p'
        write(26,*)  ws, rey_p 
      endif
      write(26,*)
c
c     transform variables to Fourier space
      call fft_r2f_2d_new(u)    !call fft_r2f_2d(u,ru)
      call fft_r2f_2d_new(v)    !call fft_r2f_2d(v,rv)
      call fft_r2f_2d_new(w)    !call fft_r2f_2d(w,rw)
      if (FLTHRM)  call fft_r2f_2d_new(tt) !call fft_r2f_2d(tt,rtt)
      if (FLTHRM2) call fft_r2f_2d_new(tt2)!call fft_r2f_2d(tt,rtt)
c
c     transform bc's (that were computed in initial)
cNOTE: to do time varying bc see subroutine advance
      call fft_r2f_2d_new_slice(ubct)
      call fft_r2f_2d_new_slice(ubcb)
      call fft_r2f_2d_new_slice(vbct)
      call fft_r2f_2d_new_slice(vbcb)
      call fft_r2f_2d_new_slice(wbct)
      call fft_r2f_2d_new_slice(wbcb)
      if (FLTHRM)  call fft_r2f_2d_new_slice(ttbct)
      if (FLTHRM)  call fft_r2f_2d_new_slice(ttbcb)
      if (FLTHRM2) call fft_r2f_2d_new_slice(tt2bct)
      if (FLTHRM2) call fft_r2f_2d_new_slice(tt2bcb)
c
      if(FLPGVAR)then
c        cfpg=2.d0*pi/ppT
        cfpg=1.d0/ppT
      else
        gradpbar=ppA
      endif
c
ccccccccccccccccccccccccccccccccccc GA ccccccccccccccccc
c      open(unit=114, file='liftx')
      open(unit=115, file='drag')
      open(unit=116, file='Sf')
      open(unit=117, file='Qf')
      open(unit=118, file='Nu')
      open(unit=119, file='Sh')
      open(unit=120, file='Ro')
      open(unit=121, file='dQfdt_A')
      open(unit=122, file='dQfdt_B')
      open(unit=125, file='dQfdt_C')
      open(unit=126, file='dQfdt_C1')
      open(unit=127, file='dQfdt_D')
      open(unit=128, file='dQfdt_DE')
      open(unit=129, file='dQfdt_E1')
      open(unit=130, file='dQfdt_E2')
      open(unit=131, file='dQfdt')
      open(unit=132, file='dSfdt')
      open(unit=133, file='dSfdt_1')
      open(unit=134, file='dSfdt_2')
      open(unit=135, file='d2Sdx2')
      open(unit=136, file='d2Sdy2')
      open(unit=137, file='Wf')
ccccccccccccccccccccccccccccccccccc GA
c     loop over the time steps
      do itime = itfirst,iend
c    
c     write the time step number to logfile
        if (mod((itime-istart),oi_cfl).eq.0) then
          write(26,*)' '
          write(26,*)' itime =  ',itime
          write(26,*)' '
        endif
c
c     perform one complete time step
        call advance(0)
        rtime  = rtime+dt
c
c     compute mean values
c        call make_means
c NOTE: adcp data for JR is output in subroutine 
c       output_adcpdata that is called in make_means above
c
c     output means
        if(FLAVER)then
c     keep these two "if" like this because oi_mean=0 -> FLAVER=.false.
c           if(mod((itime-istart),oi_mean).eq.0)call output_means
        endif
c
c     compute deposit, must be done every time step        
        if(.false..and.FLSTLE)then
           call dcopy(nx0y0,tt(1,1,nz),1,auxdepo,1) !copy bottom concentration
           call fft_f2r_2d_new_slice(auxdepo) !bring it to real space
           call comp_erosion(auxeros) !erosion returns in real space
           do j=1,ny
             do i=1,nx
                erosion(i,j)= erosion(i,j)-ws*auxeros(i,j)*dt
                deposit(i,j)= deposit(i,j)+ws*auxdepo(i,j)*dt
             enddo
           enddo
c
        endif

c     transform to real space
        if( mod((itime-istart),oi_chan).eq.0.or.
     &      mod((itime-istart),oi_cfl ).eq.0.or.
     &      mod((itime-istart),oi_mean).eq.0.or.
     &      mod((itime-istart),oi_spec).eq.0.or.
     &      mod((itime-istart),oi_gbal).eq.0.or.
     &      mod((itime-istart),oi_1d  ).eq.0)then

           call fft_f2r_2d_new(u)
           call fft_f2r_2d_new(v)
           call fft_f2r_2d_new(w)
           if (FLTHRM) call fft_f2r_2d_new(tt)
           if (FLTHRM2)call fft_f2r_2d_new(tt2)
            
c           print*,'min(u)',minval(u(1:nx,1:ny,1:nz))
c           print*,'max(u)',maxval(u(1:nx,1:ny,1:nz))
c           print*,'min(v)',minval(v(1:nx,1:ny,1:nz))
c           print*,'max(v)',maxval(v(1:nx,1:ny,1:nz))
c           print*,'min(w)',minval(w(1:nx,1:ny,1:nz))
c           print*,'max(w)',maxval(w(1:nx,1:ny,1:nz))
c           print*,'min(tt)',minval(tt(1:nx,1:ny,1:nz))
c           print*,'min(tt2)',minval(tt2(1:nx,1:ny,1:nz))
c           print*,'max(tt)',maxval(tt(1:nx,1:ny,1:nz))
c           print*,'max(tt2)',maxval(tt2(1:nx,1:ny,1:nz))
           
c NOTE: old vel are not fft, they are always in Fourier space
c
c     output a restart file
           if(mod((itime-istart),oi_chan).eq.0) call output
           if(mod((itime-istart),oi_mean).eq.0) then
             call output_2d
c             stop
           endif
c
c     output pressure (press is in Fourier space in var p allways)
           if(out_press.gt.0.and.mod((itime-istart),oi_chan).eq.0)
     &          call output_press(out_press)
c
c     compute the divergence and the cfl condition 
           if (mod((itime-istart),oi_cfl).eq.0) then
              call divg
              call courant(dt0,u,v,w)
              write(26,*)' '
c
              if(FLSTLE.and.(.not.FLGVARY))then
c     recall ws is magnitude (>0)
c     recall gz has its own sign (gz<0)
c                 igmod=sqrt(gy(1,1)**2+gz(1,1)**2)
c                 igmod=1.d0/igmod
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
                 do k=1,nz
                    do j=1,ny
                       do i=1,nx
                          rv(i,j,k)=v(i,j,k)+ws
     &                             *gy(1,k)/sqrt(gy(1,k)**2+gz(1,k)**2)
                          rw(i,j,k)=w(i,j,k)+ws
     &                             *gz(1,k)/sqrt(gy(1,k)**2+gz(1,k)**2)
                       enddo
                    enddo
                 enddo                 
C$OMP END PARALLEL DO
                 write(26,*)'  Disperse phase info:'
                 write( *,*)'  Disperse phase info:'
                 call courant(dt0,u,rv,rw)
                 write(26,*)' '
              endif
           endif
c
c     output deposition and erosion
c           if(.false..and.mod((itime-istart),oi_mean).eq.0)
c     &       call output_2d_section('deposit',deposit,erosion,deposit,2)
c
c     output Chebyshev spectrum
           if (mod((itime-istart),oi_spec).eq.0)then
              call make_chb_spec(u,espz)
              call chb_spec_out(espz,'chspu',itime)
              call make_chb_spec(v,espz)
              call chb_spec_out(espz,'chspv',itime)
              call make_chb_spec(w,espz)
              call chb_spec_out(espz,'chspw',itime)
              if(FLTHRM)then
                 call make_chb_spec(tt,espz)
                 call chb_spec_out(espz,'chspT',itime)
              endif
              if(FLTHRM2)then
                 call make_chb_spec(tt2,espz)
                 call chb_spec_out(espz,'chspT2',itime)
              endif
           endif
c
c     output 1d stats
           if (.false..and.mod((itime-istart),oi_1d).eq.0) then
              call make_output_1d(30)
              n1d = n1d + 1
           endif
c     
c     output global balances
           if (.false..and.mod((itime-istart),oi_gbal).eq.0) then
              call make_gbalrhs(KE,Ep,diss,ddtEp)
              write(31,200) rtime,KE,Ep,diss,ddtEp
              ngbal = ngbal + 1
           endif
c
c     output mean profiles for some time steps
           if(mod((itime-istart),oi_cfl).eq.0) then
              open(11,file='test-data.dat',status='unknown')
              do k=1,nz
                 write(11,11)zets(k),v(nxh,nyh,k),vm(k,1)/dble(NT(1))
     &                ,vzm(k,1)/dble(NT(1))/re
     &                ,vwm(k,1)/dble(NT(1))
     &                -vm(k,1)/dble(NT(1))*wm(k,1)/dble(NT(1))
     &                ,tt(nxh,nyh,k),ttm(k,1)/dble(NT(1))
     &                ,ttzm(k,1)/dble(NT(1))/re/pr
     &                ,wttm(k,1)/dble(NT(1))
     &                -wm(k,1)/dble(NT(1))*ttm(k,1)/dble(NT(1))
              enddo
              close(11)
           end if
c 11   FORMAT(g25.15,g25.15,g25.15,g25.15,g25.15
c     &          ,g25.15,g25.15,g25.15,g25.15,g25.15)
 11   FORMAT(10(3x,g15.8))
c
c     transform back to Fourier space
           call fft_r2f_2d_new(u)
           call fft_r2f_2d_new(v)
           call fft_r2f_2d_new(w)
           if (FLTHRM)  call fft_r2f_2d_new(tt)
           if (FLTHRM2) call fft_r2f_2d_new(tt2)
c     
        endif!outputs in real space
c
c     close and reopen gbal and 1dout when main data written, 
c     so no data is lost during a crash and restart
        if (.false..and.mod((itime-istart),oi_chan).eq.0) then
           close(30)
           open(30,file='1dout')
           rewind(30)
           do j=1,n1d
c check in make_output_1d # of variables to read
              read(30,*) diss,diss,diss,diss,diss
           enddo
           close(31)
           open(31,file='gbal')
           rewind(31)
           do j=1,ngbal
              read(31,*) diss,diss,diss,diss,diss
           enddo
        endif
c
c     if the time step is divisible by oi_timer, print timers
        if (FLTIMER.and.mod((itime-istart),oi_timer).eq.0) then
           print *, 'Intermediate timings, time = ', itime
           call print_timers
        endif
c
c     all done with one complete time step
      enddo!itime
ccccccccccccccccccccccccccccccccccccccccccccc
c      write(40,110) (zets(k),u(nxh,nyh,k),k=1,nz)
      write(40,111) (zets(k),vm( k,1)/dble(NT(1)),
     &                       ttm(k,1)/dble(NT(1)),k=1,nz)
c      write(40,111) (zets(k)+rlenz/2.d0,v(1,1,k),dsin(rtime)
c     &-exp(-(zets(k)+rlenz/2.d0)*sqrt(re/2.d0))*dsin(rtime
c     &-(zets(k)+rlenz/2.d0)*sqrt(re/2.d0)),k=1,nz)
c      write(40,111) (zets(k)+rlenz/2.d0,v(1,1,k),dsin(rtime)
c     &-exp(-(-zets(k)+rlenz/2.d0)*sqrt(re/2.d0))*dsin(rtime
c     &-(-zets(k)+rlenz/2.d0)*sqrt(re/2.d0)),k=1,nz)
ccccccccccccccccccccccccccccccccccccccccccccc
c
c     fill in and output remaining data for the final time-step
      call advance(1)
c
      close(26)
      close(30)
      close(31)
cccccccccccccccccccccccccccccccccccccccccc GA ccccccccccccccccc
c      close(unit=113)
c      close(unit=114)
      close(unit=115)
      close(unit=116)
c      close(unit=117)
c      close(unit=118)
c      close(unit=119)
c      close(unit=120)
c      close(unit=121)
c      close(unit=122)
c      close(unit=123)
c      close(unit=124)
c      close(unit=125)
ccccccccccccccccccccccccccccccccccccccccc GA
c
#ifdef FFTW3	
      call cleanfft
#endif
c
      if (FLTIMER) then
        print *, 'Final timings'
        call print_timers
      endif
c
 101  FORMAT(7i5,3f10.5)
 102  FORMAT(3f10.5,2i5)
 103  FORMAT(4i5) 
 109  FORMAT(i5,2f10.5)
 110  FORMAT(f12.8,5x,f12.8)
 111  FORMAT(f12.8,5x,f12.8,5x,f12.8)
 200  FORMAT(1p5e18.10)
c 200  FORMAT(i9,1p6e18.10)
436    END
