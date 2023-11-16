      SUBROUTINE input
#include "common.inc"
c
      INTEGER itold
      REAL*8  rtimet,dt_oldt
      INTEGER idum, i,j,k,l,n
      REAL*8  dum
      REAL*8  tt_p_2d(nx0,ny0,num_p)  ! Tempurature on all surfaces
      REAL*8  r2d(nx0,ny0,2)  ! Tempurature on all surfaces
      CHARACTER*32 anum

      open(11,file=anum(ch_file,istart),form=format_mode)
      rewind(11)
      read(11) rtime,itold,dt_old,idum,idum,idum,dum,dum,dum
      read(11) u,v,w
      read(11) pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
      read(11) pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
      close(11)
      write(26,*)'input: rtime,itold,dt_old ', rtime,itold,dt_old
      if (itold.ne.istart) then
        write(26,*)'Error in input: itold is not istart'
        stop
      endif
c
c Nadim
c Once conc.0 and conc2.0 are obtained, itime.gt.0 must be removed 
c      if (FLTHRM.and.itime.gt.0) then
      if (FLTHRM) then
        open(11,file=anum(tt_file,istart),form=format_mode)
        rewind(11)
        read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum,dum
        read(11) tt
        close(11)
        if (itold.ne.istart) then
          write(26,*)'Error in input: itold is not istart'
          stop
        endif
        if ((rtime.ne.rtimet).or.(dt_old.ne.dt_oldt)) then
          write(26,*)'Error in input: time mismatch with thermal file',
     &      rtime, rtimet, dt_old, dt_oldt
          stop
        endif
c      elseif (FLTHRM.and.itime.eq.0) then
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
c        do j=1,ny
c          do i=1,nx
c            do k=1,nz
c              tt(i,j,k) = 0.d0
cc              tt(i,j,k) = zets(k) + rlenz/2.d0 
cc              tt(i,j,k) = abs(v(i,j,k))/abs(v(1,1,1)) 
c            enddo
c          enddo
c        enddo
cC$OMP END PARALLEL DO
      endif
c
c      if (FLTHRM2.and.itime.gt.0) then
      if (FLTHRM2) then
        open(11,file=anum(tt2_file,istart),form=format_mode)
        rewind(11)
        read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum,dum
        read(11) tt2
        close(11)
        if (itold.ne.istart) then
          write(26,*)'Error in input: itold is not istart'
          stop
        endif
        if ((rtime.ne.rtimet).or.(dt_old.ne.dt_oldt)) then
          write(26,*)'Error in input: time mismatch with tt2 file',
     &      rtime, rtimet, dt_old, dt_oldt
          stop
        endif
c      elseif (FLTHRM2.and.itime.eq.0) then
cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
c        do j=1,ny
c          do i=1,nx
c            do k=1,nz
c              tt(i,j,k)  = abs(v(i,j,k))/20.0
c              tt2(i,j,k) = abs(v(i,j,k))/20.0
cc              tt2(i,j,k) = zets(k) + rlenz/2.d0
cc              tt2(i,j,k) = abs(v(i,j,k))/abs(v(1,1,1))
c            enddo
c          enddo
c        enddo
cC$OMP END PARALLEL DO
      endif
c
      if (FLSURF) then
        open(11,file=anum(lag_file,istart),form=format_mode)
        rewind(11)
        read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum,dum
        read(11) z2d,r2d,tt_p_2d 
        close(11)
        if (itold.ne.istart) then
          write(26,*)'Error in input: itold is not istart'
          stop
        endif
        if ((rtime.ne.rtimet).or.(dt_old.ne.dt_oldt)) then
          write(26,*)'Error in input: time mismatch with lag file',
     &      rtime, rtimet, dt_old, dt_oldt
          stop
        endif

c        z_o_avg = sum(z2d(1:nx,1:ny))/dble(nx*ny)+1.d0-rlenz/2.d0 

cC$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
c        do j=1,ny
c          do i=1,nx
c            z2d(i,j)= z2d(i,j) - z_o_avg 
c          enddo
c        enddo
cC$OMP END PARALLEL DO

        do n=1,3,2 ! Nadim Must verify nx*ny > nxyp*nzp
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
          do i=1,nx
            tt_p((i-1)*ny+j,n) = tt_p_2d(i,j,n)
          enddo
        enddo
C$OMP END PARALLEL DO
        enddo

        do n=5,7,2 ! Nadim Must verify nx*ny > nxyp*nzp
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k,l)
        do k=1,nzp
          do l=1,nxyp
            tt_p((l-1)*nzp+k,n) = tt_p_2d(l,k,n) 
          enddo
        enddo
C$OMP END PARALLEL DO
        enddo

C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k,l)
        do k=1,nzp
          do l=1,nxyp
            rf( (l-1)*nzp+k,1) = r2d(l,k,1) 
            rf( (l-1)*nzp+k,2) = r2d(l,k,2) 
          enddo
        enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
          do i=1,nx
            z_o( (i-1)*ny+j,1) = z2d(    i,j,1) 
            z_o( (i-1)*ny+j,3) = z2d(    i,j,2) 
          enddo
        enddo
C$OMP END PARALLEL DO
c        
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k,l)
        do k=1,nzp
          do l=1,nxyp
            x_o( (l-1)*nzp+k,5) = rf((l-1)*nzp+k,1) 
     &              *dcos(dble(l-1)*deltaq)+x01
            if(x_o( (l-1)*nzp+k,5).lt.0) then
                    x_o( (l-1)*nzp+k,5) = x_o( (l-1)*nzp+k,5) + rlenx
            elseif(x_o( (l-1)*nzp+k,5).gt.rlenx) then
                    x_o( (l-1)*nzp+k,5) = x_o( (l-1)*nzp+k,5) - rlenx
            endif

            y_o( (l-1)*nzp+k,5) = rf((l-1)*nzp+k,1) 
     &              *dsin(dble(l-1)*deltaq)+y01
            if(y_o( (l-1)*nzp+k,5).lt.0) then
                    y_o( (l-1)*nzp+k,5) = y_o( (l-1)*nzp+k,5) + rleny
            elseif(y_o( (l-1)*nzp+k,5).gt.rleny) then
                    y_o( (l-1)*nzp+k,5) = y_o( (l-1)*nzp+k,5) - rleny
            endif

            x_o( (l-1)*nzp+k,7) = rf((l-1)*nzp+k,2) 
     &              *dcos(dble(l-1)*deltaq)+x02
            if(x_o( (l-1)*nzp+k,7).lt.0) then
                    x_o( (l-1)*nzp+k,7) = x_o( (l-1)*nzp+k,7) + rlenx
            elseif(x_o( (l-1)*nzp+k,7).gt.rlenx) then
                    x_o( (l-1)*nzp+k,7) = x_o( (l-1)*nzp+k,7) - rlenx
            endif

            y_o( (l-1)*nzp+k,7) = rf((l-1)*nzp+k,2) 
     &              *dsin(dble(l-1)*deltaq)+y02
            if(y_o( (l-1)*nzp+k,7).lt.0) then
                    y_o( (l-1)*nzp+k,7) = y_o( (l-1)*nzp+k,7) + rleny
            elseif(y_o( (l-1)*nzp+k,7).gt.rleny) then
                    y_o( (l-1)*nzp+k,7) = y_o( (l-1)*nzp+k,7) - rleny
            endif
          enddo
        enddo
C$OMP END PARALLEL DO
      else !FLSURF
        do n=1,num_p,2 ! Nadim Must verify nx*ny > nxyp*nzp
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
          do i=1,nx
c            z_o( (i-1)*ny+j,1) = rlenz/2.d0 - 1.d0 
            tt_p((i-1)*ny+j,n) = 0.d0 
          enddo
        enddo
C$OMP END PARALLEL DO
        enddo
      endif !FLSURF
c
c if ws.ne.0 input deposit
      if(.false..and.FLSTLE)then 
         open(11,file=anum('deposit',istart),form=format_mode)
         rewind(11)
         read(11) rtimet,itold,dt_oldt,idum,idum,idum,dum,dum,dum,idum
         read(11) deposit,erosion
         close(11)
      endif
c
      if(iters.lt.0) then
        call output
        stop
      endif

      return
      END

      SUBROUTINE output
#include "common.inc"
c
      CHARACTER*32 anum
      INTEGER      n,i,j
      REAL*8  tt_p_2d(nx0,ny0,num_p) 
      REAL*8  r2d(nx0,ny0,2) 
      REAL*8 dttpdn_2d(nx0,ny0,num_p)
      REAL*8 dtt2pdn_2d(nx0,ny0,num_p)
c
      if (FLTIMER) call begin_timer(30)
      if(iters.lt.0) then
        open(21,file=anum(ch_file,0),form=format_mode)
        write(21) 0.d0,0,dt,nx0,ny0,nz,deltax,deltay,rlenz
      else
        open(21,file=anum(ch_file,itime),form=format_mode)
        write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      endif
      write(21) u,v,w
      write(21) pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
      write(21) pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
      close(21)
c
      if(FLMOV)then
        do n=1,num_p,2
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
          do i=1,nx
            tt_p_2d(i,j,n) = tt_p((i-1)*ny+j,n) 
            dttpdn_2d(i,j,n)  = dttp_tildedn( (i-1)*ny+j,n)
            dtt2pdn_2d(i,j,n) = dtt2p_tildedn((i-1)*ny+j,n)
          enddo
        enddo
C$OMP END PARALLEL DO
        enddo
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
          do i=1,nx
            z2d(i,j,1) = z_o((i-1)*ny+j,1)
            z2d(i,j,2) = z_o((i-1)*ny+j,3)
          enddo
        enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j=1,nzp
        do i=1,nxyp
          r2d(i,j,1) = rf((i-1)*nzp+j,1)
          r2d(i,j,2) = rf((i-1)*nzp+j,2)
        enddo
      enddo
C$OMP END PARALLEL DO


        if(iters.lt.0) then
          open(21,file=anum(lag_file,0),form=format_mode)
          write(21) 0.d0,0,dt,nx0,ny0,nz,deltax,deltay,rlenz
        else
          open(21,file=anum(lag_file,itime),form=format_mode)
          write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
        endif
        write(21) z2d, r2d, tt_p_2d, dttpdn_2d, dtt2pdn_2d
        close(21)
      endif
c     
      if(FLTHRM)then
      if(iters.lt.0) then
        open(21,file=anum(tt_file,0),form=format_mode)
        write(21) 0.d0,0,dt,nx0,ny0,nz,deltax,deltay,rlenz
      else
        open(21,file=anum(tt_file,itime),form=format_mode)
        write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      endif
        write(21) tt
        close(21)
      endif
c     
      if(FLTHRM2)then
      if(iters.lt.0) then
        open(21,file=anum(tt2_file,0),form=format_mode)
        write(21) 0.d0,0,dt,nx0,ny0,nz,deltax,deltay,rlenz
      else
        open(21,file=anum(tt2_file,itime),form=format_mode)
        write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      endif
        write(21) tt2
        close(21)
      endif
c     
      if (FLTIMER) call end_timer(30)
      return
      END
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                          Nadim Zgheib                            C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      SUBROUTINE output_2d
#include "common.inc"
      REAL*8 val2d1(ny0,nz),val2d2(ny0,nz),val2d3(ny0,nz),
     &          val2dtt(ny0,nz), val2dtt2(ny0,nz)
      INTEGER i,j,k,n,l
      REAL*8 tt_p_2d(nx0,ny0,num_p),r2d(nx0,ny0,2)
      REAL*8 dttpdn_2d(nx0,ny0,num_p)
      REAL*8 dtt2pdn_2d(nx0,ny0,num_p)

      CHARACTER*32 anum

      if (FLTIMER) call begin_timer(36)

        do n=1,3,2
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
          do i=1,nx
            tt_p_2d(i,j,n)    = tt_p(         (i-1)*ny+j,n)
            dttpdn_2d(i,j,n)  = dttp_tildedn( (i-1)*ny+j,n)
            dtt2pdn_2d(i,j,n) = dtt2p_tildedn((i-1)*ny+j,n)
          enddo
        enddo
C$OMP END PARALLEL DO
        enddo
c
        do n=5,7,2
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j=1,nzp
        do i=1,nxyp
            tt_p_2d(i,j,n)    = tt_p(         (i-1)*nzp+j,n) 
            dttpdn_2d(i,j,n)  = dttp_tildedn( (i-1)*nzp+j,n) 
            dtt2pdn_2d(i,j,n) = dtt2p_tildedn((i-1)*nzp+j,n) 
          enddo
        enddo
C$OMP END PARALLEL DO
        enddo
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
        do j=1,ny
          do i=1,nx
            z2d(i,j,1) = z_o((i-1)*ny+j,1)
            z2d(i,j,2) = z_o((i-1)*ny+j,3)
          enddo
        enddo
C$OMP END PARALLEL DO
c
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j)
      do j=1,nzp
        do i=1,nxyp
          r2d(i,j,1) = rf((i-1)*nzp+j,1)
          r2d(i,j,2) = rf((i-1)*nzp+j,2)
        enddo
      enddo
C$OMP END PARALLEL DO

      if((FLMOV.and.itime.eq.0).or.
     &        (FLMOV.and.mod((itime-istart),oi_chan).ne.0))then
        open(21,file=anum(lag_file,itime),form=format_mode)
        write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
        write(21) z2d, r2d, tt_p_2d, dttpdn_2d, dtt2pdn_2d
        close(21)
      endif
c
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,k)
      do k=1,nz
        do j=1,ny
          val2d1(j,k) = u(nx/2,j,k)
          val2d2(j,k) = v(nx/2,j,k)
          val2d3(j,k) = w(nx/2,j,k)
        enddo
      enddo
C$OMP END PARALLEL DO

      if(FLTHRM)then
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,k)
         do k=1,nz
           do j=1,ny
             val2dtt(j,k) = tt(nx/2,j,k)
           enddo
         enddo
C$OMP END PARALLEL DO
      endif

      if(FLTHRM2)then
C$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(j,k)
         do k=1,nz
           do j=1,ny
             val2dtt2(j,k) = tt2(nx/2,j,k)
           enddo
         enddo
C$OMP END PARALLEL DO
      endif
cc
cc write averaged value
c      open(21,file=anum('vel_xy',itime),form=format_mode)
c      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
c      write(21) valxy1,valxy2,valxy3
c      close(21)
cc
c write averaged value
      if(FLTHRM.and.FLTHRM2)then
        if(iters.lt.0) then
          open(21,file=anum('vel_tt_tt2_yz',0),form=format_mode)
          write(21) 0.d0,0,dt,nx0,ny0,nz,deltax,deltay,rlenz
        else
          open(21,file=anum('vel_tt_tt2_yz',itime),form=format_mode)
          write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
        endif
        write(21) val2d1,val2d2,val2d3,val2dtt,val2dtt2
        close(21)
      elseif(FLTHRM)then
        if(iters.lt.0) then
          open(21,file=anum('vel_tt_yz',0),form=format_mode)
          write(21) 0.d0,0,dt,nx0,ny0,nz,deltax,deltay,rlenz
        else
          open(21,file=anum('vel_tt_yz',itime),form=format_mode)
          write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
        endif
        write(21) val2d1,val2d2,val2d3,val2dtt
        close(21)
      else
        if(iters.lt.0) then
          open(21,file=anum('vel_yz',0),form=format_mode)
          write(21) 0.d0,0,dt,nx0,ny0,nz,deltax,deltay,rlenz
        else
          open(21,file=anum('vel_yz',itime),form=format_mode)
          write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
        endif
        write(21) val2d1,val2d2,val2d3
        close(21)
      endif

      if (FLTIMER) call end_timer(36)
      return
      END

*=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      SUBROUTINE output_2d_section(name,val2d_1,val2d_2,val2d_3,dim)
#include "common.inc"
      REAL*8 val2d_1(nx0,ny0),val2d_2(nx0,ny0),val2d_3(nx0,ny0)
      INTEGER dim
      CHARACTER*(*) name
c
      CHARACTER*32 anum
c
      if (FLTIMER) call begin_timer(36)
c write averaged value
      open(21,file=anum(name,itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz,0
      if(dim.eq.1) write(21) val2d_1
      if(dim.eq.2) write(21) val2d_1,val2d_2
      if(dim.eq.3) write(21) val2d_1,val2d_2,val2d_3
      close(21)
c
      if (FLTIMER) call end_timer(36)
      return
      END
c
      SUBROUTINE input_means
#include "common.inc"
c
      CHARACTER*32 anum
      INTEGER nz_old,np_old,idum
      REAL*8 rlenzaux,dum
c
c      if (FLTIMER) call begin_timer()
      open(21,file=anum('means',itime),form=format_mode)
      read(21)dum,idum,dum,idum,idum,nz_old,np_old,dum,dum,rlenzaux,NT
      read(21)dum,dum,dum
      if(nz.ne.nz_old.or.np.ne.np_old.or.rlenz.ne.rlenzaux)then
         write(26,*)'Error in input of mean: nz_old ne nz'
         write(26,*)'                        or rlenz_old ne rlenz'
         write(26,*)'                        or np_old ne np'
         stop
      endif
      read(21)um,vm,wm,ttm,u2m,v2m,w2m,tt2m,u3m,v3m,w3m,tt3m,u4m,v4m
     &,w4m,tt4m,uvm,uwm,uttm,vwm,vttm,wttm,uxm,uym,uzm,ux2m,uy2m,uz2m
     &,ux3m,uy3m,uz3m,ux4m,uy4m,uz4m,vxm,vym,vzm,vx2m,vy2m,vz2m,vx3m
     &,vy3m,vz3m,vx4m,vy4m,vz4m,wxm,wym,wzm,wx2m,wy2m,wz2m,wx3m,wy3m
     &,wz3m,wx4m,wy4m,wz4m,ttxm,ttym,ttzm,ttx2m,tty2m,ttz2m,ttx3m,tty3m
     &,ttz3m,ttx4m,tty4m,ttz4m
      read(21)uespxm,uespym,vespxm,vespym,wespxm,wespym
      read(21)pwm,pm,u2wm,v2wm,uuzm,uwxm,vvzm,vwym,wwzm
      close(21)
c
c      if (FLTIMER) call end_timer()
      return
      END
c
      SUBROUTINE output_means
#include "common.inc"
c
      CHARACTER*32 anum
c
c      if (FLTIMER) call begin_timer()
      open(21,file=anum('means',itime),form=format_mode)
      write(21)rtime,itime,dt,nx0,ny0,nz,np,deltax,deltay,rlenz,NT
      write(21)re,pr,ra
      write(21)um,vm,wm,ttm,u2m,v2m,w2m,tt2m,u3m,v3m,w3m,tt3m,u4m,v4m
     &,w4m,tt4m,uvm,uwm,uttm,vwm,vttm,wttm,uxm,uym,uzm,ux2m,uy2m,uz2m
     &,ux3m,uy3m,uz3m,ux4m,uy4m,uz4m,vxm,vym,vzm,vx2m,vy2m,vz2m,vx3m
     &,vy3m,vz3m,vx4m,vy4m,vz4m,wxm,wym,wzm,wx2m,wy2m,wz2m,wx3m,wy3m
     &,wz3m,wx4m,wy4m,wz4m,ttxm,ttym,ttzm,ttx2m,tty2m,ttz2m,ttx3m,tty3m
     &,ttz3m,ttx4m,tty4m,ttz4m
      write(21)uespxm,uespym,vespxm,vespym,wespxm,wespym
      write(21)pwm,pm,u2wm,v2wm,uuzm,uwxm,vvzm,vwym,wwzm
      close(21)
c
c      if (FLTIMER) call end_timer()
      return
      END
c
      SUBROUTINE mean5_out(cl1,cl2,cl3,cl4,cl5,m1,m2,m3,m4,m5,outfile,t)
#include "common.inc"
      INTEGER t
      REAL*8 cl1(nz),cl2(nz),cl3(nz),cl4(nz),cl5(nz),m1,m2,m3,m4,m5
      CHARACTER*(*) outfile
c
      INTEGER k
      CHARACTER*32 anum
c
      open(21,file=anum(outfile,t))
      write(21,1000) m1, m2, m3, m4, m5
      do k=1,nz
        write(21,2000) zets(k), cl1(k), cl2(k), cl3(k), cl4(k), cl5(k)
      enddo
      close(21)
c
      return
 1000 FORMAT('#   means:        ', 1p9e18.10)
 2000 FORMAT(1p9e18.10)
      END
c
      SUBROUTINE spec_out(espx,dspx,espy,dspy,outfile,t)
#include "common.inc"
      INTEGER t
      REAL*8 espx(nxh,nz),dspx(nxh,nz),espy(nyh,nz),dspy(nyh,nz)
      CHARACTER*(*) outfile
c
      CHARACTER*32 anum
c
      open(21,file=anum(outfile,t),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      write(21) espx
      write(21) dspx
      write(21) espy
      write(21) dspy
      close(21)
c
      return
      END
c
      SUBROUTINE chb_spec_out(espz,outfile,t)
#include "common.inc"
      INTEGER t
      REAL*8 espz(nz)
      CHARACTER*(*) outfile
c
      CHARACTER*32 anum
c
      open(21,file=anum(outfile,t),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      write(21) espz
      close(21)
c
      return
      END
c
      SUBROUTINE output_press(code)
#include "common.inc"
      INTEGER code
c
c     This subroutine postprocesses and outputs pressure for output 
c     according to code:
c     code=1 -> output dyn press + hydro press due to strat by variable tt
c     code=2 -> output dyn press
c     Needs tt in physical space
c
      real*8 tmeanp,hydp,s
      integer i,j,k,l
      CHARACTER*32 anum
c
c      if (FLTIMER) call begin_timer(36)
c
c     copy p to rw
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(k)
      do k=1,nz
        call dcopy(nx0y,p(1,1,k),1,rw(1,1,k),1)
      enddo
C$OMP END PARALLEL DO
c
c     IFFT rw -> rw is pressure
      call fft_f2r_2d_new(rw)
c
c     substract stratification pressure 
c     and compute total mean dynamic pressure
      tmeanp = 0.d0
      if(code.eq.2)then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,l,hydp,s)
C$OMP&                           ,REDUCTION(+:tmeanp)
         do k=1,nz
            s=0.d0
            do j=1,ny
               do i=1,nx
                  hydp=0.d0
                  do l=1,k
                     hydp=hydp+tt(i,j,l)*(rlenz*wzet(l))
                  enddo
                  rw(i,j,k)=rw(i,j,k)-hydp
                  s=s+rw(i,j,k)
               enddo
            enddo
            s=s/dble(nx*ny)
            tmeanp=tmeanp+s*wzet(k)
         enddo
C$OMP END PARALLEL DO
      elseif(code.eq.1)then
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k,s)
C$OMP&                           ,REDUCTION(+:tmeanp)
         do k=1,nz
            s=0.d0
            do j=1,ny
               do i=1,nx
                  s=s+rw(i,j,k)
               enddo
            enddo
            s=s/dble(nx*ny)
            tmeanp=tmeanp+s*wzet(k)
         enddo
C$OMP END PARALLEL DO
      endif
c
c     substract mean pressure to normalize output
C$OMP PARALLEL DO DEFAULT(SHARED),PRIVATE(i,j,k)
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  rw(i,j,k)=rw(i,j,k)-tmeanp
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO
c
      open(21,file=anum('press',itime),form=format_mode)
      write(21) rtime,itime,dt,nx0,ny0,nz,deltax,deltay,rlenz
      write(21) rw
      close(21)
c
c      if (FLTIMER) call end_timer(36)
      return
      END
c
      SUBROUTINE output_adcpdata
C     NOTE: this subroutine assumes physical space data in ru,rv,rw and rtt
#include "common.inc"
      integer nxj,nyj
      parameter(nxj=21,nyj=21)
      real   ucil(nxj,nyj,nz),vcil(nxj,nyj,nz),wcil(nxj,nyj,nz)
      real   uslicel(3,ny,nz),vslicel(3,ny,nz),wslicel(3,ny,nz)
      real   uslicet(nx,3,nz),vslicet(nx,3,nz),wslicet(nx,3,nz)
      integer i,j,k
      CHARACTER*32 anum
ccccccccccccccccccccccccc
c      real   ucm(nz),vcm(nz),wcm(nz)
ccccccccccccccccccccccccc
c
      do k=1,nz
         do j=1,nyj
            do i=1,nxj
               ucil(i,j,k)=ru(nxhp-nxj/2-1+i,nyhp-nyj/2-1+j,k)
               vcil(i,j,k)=rv(nxhp-nxj/2-1+i,nyhp-nyj/2-1+j,k)
               wcil(i,j,k)=rw(nxhp-nxj/2-1+i,nyhp-nyj/2-1+j,k)
            enddo
         enddo
      enddo
ccccccccccccccccccccccccc
c      do k=1,nz
c         ucm(k)=0.0
c         vcm(k)=0.0
c         wcm(k)=0.0
c      enddo
c      do k=1,nz
c         do j=1,nyj
c            do i=1,nxj
c               ucm(k)=ucm(k)+ucil(i,j,k)/dble(nxj*nyj)
c               vcm(k)=vcm(k)+vcil(i,j,k)/dble(nxj*nyj)
c               wcm(k)=wcm(k)+wcil(i,j,k)/dble(nxj*nyj)
c            enddo
c         enddo
c      enddo
c      write(12,*)(zets(k),ucm(k),vcm(k),wcm(k),z=1,nz)
c      write(12,*)(zets(k),ucil(1,1,k),vcil(1,1,k),wcil(1,1,k),k=1,nz)
ccccccccccccccccccccccccc
c
      do k=1,nz
         do j=1,ny
            do i=1,3
               uslicel(i,j,k)=ru(nxhp-2+i,j,k)
               vslicel(i,j,k)=rv(nxhp-2+i,j,k)
               wslicel(i,j,k)=rw(nxhp-2+i,j,k)
            enddo
         enddo
      enddo
c
      do k=1,nz
         do j=1,3
            do i=1,nx
               uslicet(i,j,k)=ru(i,nyhp-2+j,k)
               vslicet(i,j,k)=rv(i,nyhp-2+j,k)
               wslicet(i,j,k)=rw(i,nyhp-2+j,k)
            enddo
         enddo
      enddo
c
      open(21,file=anum('vel-cil',itime),form=format_mode)
      write(21) nxj,nyj,nz,re,rlenx,rleny,rlenz
      write(21) ucil,vcil,wcil
      close(21)
c
      open(21,file=anum('vel-slicel',itime),form=format_mode)
      write(21) 3,ny,nz,re,rlenx,rleny,rlenz
      write(21) uslicel,vslicel,wslicel
      close(21)
c      
      open(21,file=anum('vel-slicet',itime),form=format_mode)
      write(21) nx,3,nz,re,rlenx,rleny,rlenz
      write(21) uslicet,vslicet,wslicet
      close(21)
c      
      return
      END      
c
      FUNCTION anum(string,num)
      IMPLICIT none
      CHARACTER*32 anum
      CHARACTER*(*) string
      INTEGER num
c
      INTEGER i,n,dot,fin
c
      anum = string
      dot = index(anum,' ')
      if (dot.gt.0.and.num.ge.0) then
        if (dot.le.32) anum(dot:dot) = '.'
        fin = dot
        n = num
 100    fin = fin + 1
        n = n/10
        if (n.gt.0) go to 100
        n = num
        do i=fin,dot+1,-1
          if (i.le.32) anum(i:i) = char(ichar('0')+mod(n,10))
          n = n/10
        enddo
      endif
c
      return
      END
