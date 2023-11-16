      IMPLICIT none
#include "global.inc"
#include "fft.inc"
#include "flags.inc"
#include "timers.inc"
c
      INTEGER istart,itime,iters
      INTEGER oi_chan,out_press
      INTEGER oi_timer,oi_cfl,oi_spec,oi_mean,oi_gbal,oi_1d
      INTEGER msvx, msvy
      REAL*8 re,rlenx,rleny,rlenz,deltax,deltay,idx,idy
      REAL*8 dt,dt_old,rtime,pi,ipin,cflmax,gradpbar,ppA,ppT,cfpg
     &       ,theta
      REAL*8 pr,ra,Tmax 
      REAL*8 uvab,   uvat,   uvbb,   uvbt,    ugb,    ugt,vgb,vgt
      REAL*8  wab,    wat,    wbb,    wbt,    wgb,    wgt
      REAL*8 ttab(m),ttat(m),ttbb(m),ttbt(m),ttgb(m),ttgt(m)
      REAL*8 ws,rey_p
      REAL*8 cferos1,cferos2 
      REAL*8 cfnl1v(3),cfnl2v(3),cfdifv(3),cfdelv(3),cfcumv(3)
      REAL*8 zets(nz),deltaz(nzm),idz(nzm),ndeltaz(nzm),nidz(nzm),
     &  wzet(nz)
      REAL*8 wavex(nx),wavexs(nx),wavey(ny),waveys(ny)
      REAL*8 chb(nz0,nz)
      REAL*8 chbd1(nz0,nz),chbd2(nz0,nz)
      REAL*8 cn_drdr(nz0,nz),ci_drdr(nz0,nz),e_drdr(nz)
      REAL*8 uvcn_mimi(nzmm0,nzmm),uvci_mimi(nzmm0,nzmm),uve_mimi(nzmm)
      REAL*8 wcn_mimi(nzmm0,nzmm),wci_mimi(nzmm0,nzmm),we_mimi(nzmm)
      REAL*8 ttcn_mimi(nzmm0,nzmm),ttci_mimi(nzmm0,nzmm),tte_mimi(nzmm)
      REAL*8 cn_nodr(nz0,nz),ci_nodr(nz0,nz),e_nodr(nz)
      REAL*8 cn_nono(nzmm0,nzmm),ci_nono(nzmm0,nzmm),e_nono(nzmm)
      REAL*8 mlt(nz0,nx0,nyh)
      REAL*8 pbxnb(nx0,nyh),pbynb(nx0,nyh)
      REAL*8 pbxnt(nx0,nyh),pbynt(nx0,nyh)
      REAL*8 pbxob(nx0,nyh),pbyob(nx0,nyh)
      REAL*8 pbxot(nx0,nyh),pbyot(nx0,nyh)
      REAL*8 pcxnb(nx0,nyh),pcynb(nx0,nyh)
      REAL*8 pcxnt(nx0,nyh),pcynt(nx0,nyh)
      REAL*8 pcxob(nx0,nyh),pcyob(nx0,nyh)
      REAL*8 pcxot(nx0,nyh),pcyot(nx0,nyh)
      REAL*8   prbc1(nzmm),  prbcn(nzmm),
     &       uvprbc1(nzmm), uvprbcn(nzmm),
     &        wprbc1(nzmm), wprbcn(nzmm),
     &       ttprbc1(nzmm),ttprbcn(nzmm)
      REAL*8 u(nx0,ny0,nz),v(nx0,ny0,nz),w(nx0,ny0,nz),p(nx0,ny0,nz)
      REAL*8 ru(nx0,ny0,nz),rv(nx0,ny0,nz),rw(nx0,ny0,nz),rp(nx0,ny0,nz)
      REAL*8 h1(nx0,ny0,nz),h2(nx0,ny0,nz),h3(nx0,ny0,nz)
c      REAL*8 g(nx0,ny0,nz)
      REAL*8 auxt1(nz0,nx0,nyh),auxt2(nz0,nx0,nyh)
      REAL*8 tt(nx0,ny0,nz),rtt(nx0,ny0,nz),htt(nx0,ny0,nz)
      REAL*8 u_yz_plane(ny0,2),v_yz_plane(ny0,2),w_yz_plane(ny0,2)
      REAL*8 u_2d_yz_plane(ny0,nz),v_2d_yz_plane(ny0,nz)
     &      ,w_2d_yz_plane(ny0,nz)
      REAL*8 ul(nxl0,nyl0,nz),vl(nxl0,nyl0,nz),wl(nxl0,nyl0,nz)
      REAL*8 ddxul(nxl0,nyl0,nz),ddyul(nxl0,nyl0,nz),ddzul(nxl0,nyl0,nz)
      REAL*8 svx(nx),svy(nyh)
      REAL*8 epssvx,epssvy
      REAL*8 deposit(nx0,ny0),erosion(nx0,ny0)
      REAL*8 ubct(nx0,ny0),ubcb(nx0,ny0),
     &       ubctax(nx0,ny0),ubcbax(nx0,ny0)
      REAL*8 vbct(nx0,ny0),vbcb(nx0,ny0),
     &       vbctax(nx0,ny0),vbcbax(nx0,ny0),
     &       uvmbct1,uvmbcb1,uvmbctn,uvmbcbn
      REAL*8 wbct(nx0,ny0),wbcb(nx0,ny0),
     &       wmbct1,wmbcb1,wmbctn,wmbcbn
      REAL*8 ttbct(nx0,ny0),ttbcb(nx0,ny0),
     &       ttmbct1,ttmbcb1,ttmbctn,ttmbcbn
      REAL*8 gy(nyl0,nz),gz(nyl0,nz)
      CHARACTER*32 ch_file,tt_file,format_mode,lag_file
c
      INTEGER kts(11)
      INTEGER NT(np)
      REAL*8  um(nz,np),vm(nz,np),wm(nz,np),ttm(nz,np)
     &       ,u2m(nz,np),v2m(nz,np),w2m(nz,np),tt2m(nz,np)
     &       ,u3m(nz,np),v3m(nz,np),w3m(nz,np),tt3m(nz,np)
     &       ,u4m(nz,np),v4m(nz,np),w4m(nz,np),tt4m(nz,np)
     &       ,uvm(nz,np),uwm(nz,np),uttm(nz,np)
     &       ,vwm(nz,np),vttm(nz,np)
     &       ,wttm(nz,np)
      REAL*8  uxm( nz,np),uym( nz,np),uzm( nz,np)
     &       ,ux2m(nz,np),uy2m(nz,np),uz2m(nz,np)
     &       ,ux3m(nz,np),uy3m(nz,np),uz3m(nz,np)
     &       ,ux4m(nz,np),uy4m(nz,np),uz4m(nz,np)
      REAL*8  vxm( nz,np),vym( nz,np),vzm( nz,np)
     &       ,vx2m(nz,np),vy2m(nz,np),vz2m(nz,np)
     &       ,vx3m(nz,np),vy3m(nz,np),vz3m(nz,np)
     &       ,vx4m(nz,np),vy4m(nz,np),vz4m(nz,np)
      REAL*8  wxm( nz,np),wym( nz,np),wzm( nz,np)
     &       ,wx2m(nz,np),wy2m(nz,np),wz2m(nz,np)
     &       ,wx3m(nz,np),wy3m(nz,np),wz3m(nz,np)
     &       ,wx4m(nz,np),wy4m(nz,np),wz4m(nz,np)
      REAL*8  ttxm( nz,np),ttym( nz,np),ttzm( nz,np)
     &       ,ttx2m(nz,np),tty2m(nz,np),ttz2m(nz,np)
     &       ,ttx3m(nz,np),tty3m(nz,np),ttz3m(nz,np)
     &       ,ttx4m(nz,np),tty4m(nz,np),ttz4m(nz,np)
      REAL*8  uespxm(nxh,nz,np),uespym(nyh,nz,np)  
ccc     &       ,udspxm(nxh,nz,np),udspym(nyh,nz,np)
      REAL*8  vespxm(nxh,nz,np),vespym(nyh,nz,np)  
ccc     &       ,vdspxm(nxh,nz,np),vdspym(nyh,nz,np)
      REAL*8  wespxm(nxh,nz,np),wespym(nyh,nz,np)  
ccc     &       ,wdspxm(nxh,nz,np),wdspym(nyh,nz,np)
      REAL*8  pwm(nz,np),pm(nz,np)
     &       ,u2wm(nz,np),v2wm(nz,np)
     &       ,uuzm(nz,np),uwxm(nz,np)
     &       ,vvzm(nz,np),vwym(nz,np)
     &       ,wwzm(nz,np)
c
ccccccccccccccccccccccc 2nd scalar field ccccccccccccccccccccccc
      REAL*8 pr2,ra2,Tmax2
      REAL*8 epssvx2,epssvy2
      INTEGER msvx2, msvy2
      REAL*8 ul2(nxl0,nyl0,nz),vl2(nxl0,nyl0,nz),wl2(nxl0,nyl0,nz)
      REAL*8 svx2(nx),svy2(nyh)
      REAL*8 ru2(nx0,ny0,nz),rv2(nx0,ny0,nz),rw2(nx0,ny0,nz),
     &       rp2(nx0,ny0,nz)
      REAL*8 ddxul2(nxl0,nyl0,nz)

      REAL*8 tt2ab(m),tt2at(m),tt2bb(m),tt2bt(m),tt2gb(m),tt2gt(m)
      REAL*8 tt2cn_mimi(nzmm0,nzmm),tt2ci_mimi(nzmm0,nzmm)
     &      ,tt2e_mimi(nzmm)
      REAL*8 tt2prbc1(nzmm),tt2prbcn(nzmm)
      REAL*8 tt2(nx0,ny0,nz),rtt2(nx0,ny0,nz),htt2(nx0,ny0,nz)
      REAL*8 tt2bct(nx0,ny0),tt2bcb(nx0,ny0),
     &       tt2mbct1,tt2mbcb1,tt2mbctn,tt2mbcbn
      CHARACTER*32 tt2_file
      REAL*8  forcing_tt2(nx0,ny0,nz)
      REAL*8  dtt2_c_dum_n(n_ll,num_p),tt2p_tilde_n(n_ll,num_p)
ccccccccccccccccccccccc 2nd scalar field ccccccccccccccccccccccc
c

ccccccccccccccccccccccc Immersed boundary cccccccccccccccccccccc
      REAL*8  tt_p(n_ll,num_p),tt2_p(n_ll,num_p)
      REAL*8  dttp_tildedn(n_ll,num_p),dtt2p_tildedn(n_ll,num_p)
      REAL*8  T_wall, Y_wall, T_inf, Y_inf 
      REAL*8  z_o_avg, AF, diff_coef
      INTEGER n_l(num_p), num_x_o
      REAL*8  deltaq, Rpillar
      REAL*8  x01, y01, x02, y02
ccccccccccccccccccccccc Immersed boundary cccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccc GA cccccccccccccccccc
c     common variables used by the ibm
      REAL*8  xets(nx0),yets(ny0)
      REAL*8  zplus,rplus

      REAL*8  z2d(nx0,ny0,2)
      REAL*8  rf(n_ll,2)
      REAL*8  exner_rhs(n_ll,num_p)
      REAL*8  forcing_x(nx0,ny0,nz),forcing_y(nx0,ny0,nz)
     &       ,forcing_z(nx0,ny0,nz),forcing_tt(nx0,ny0,nz)
      REAL*8  ddf(n_ll,num_p,4,4,6)
      INTEGER irkk
      REAL*8  x_o(n_ll,num_p),y_o(n_ll,num_p),z_o(n_ll,num_p)
     &       ,dv_l(n_ll,num_p)
      REAL*8  z_o_nm1(n_ll), rho_fro_nm1(n_ll)
      REAL*8  p_iw(n_ll,num_p),p_ie(n_ll,num_p)
     &       ,p_js(n_ll,num_p),p_jn(n_ll,num_p)
     &       ,p_kb(n_ll,num_p),p_kt(n_ll,num_p)
      REAL*8  p_iw_n(n_ll,num_p),p_ie_n(n_ll,num_p)
     &       ,p_js_n(n_ll,num_p),p_jn_n(n_ll,num_p)
     &       ,p_kb_n(n_ll,num_p),p_kt_n(n_ll,num_p)
      REAL*8  u_p(n_ll,num_p),v_p(n_ll,num_p),w_p(n_ll,num_p)
      REAL*8  chbd1_ibm(nz,nz)
      REAL*8  dpx(nx0,ny0,nz),dpy(nx0,ny0,nz),dpz(nx0,ny0,nz),drag,lift
cccccccccccccccccccccccccccccccccccccccccc GA
c
      COMMON /iparm/   istart,itime,iters
      COMMON /oiparm/  oi_chan,oi_timer,out_press
     &                ,oi_cfl,oi_spec,oi_mean,oi_gbal,oi_1d
      COMMON /rparm1/  re,rlenx,rleny,rlenz,deltax,deltay,idx,idy
      COMMON /rparm2/  dt,dt_old,rtime,pi,ipin,cflmax,gradpbar
     &                                  ,ppA,ppT,cfpg,theta
      COMMON /rparmt/  pr,ra,Tmax 
      COMMON /abg/     uvab,uvat,uvbb,uvbt, ugb, ugt,
     &                                      vgb, vgt,
     &                  wab, wat, wbb, wbt, wgb, wgt,
     &                 ttab,ttat,ttbb,ttbt,ttgb,ttgt
      COMMON /rparmw/  ws,rey_p
     &                 ,cferos1,cferos2 
      COMMON /rkcv/    cfnl1v,cfnl2v,cfdifv,cfdelv,cfcumv
      COMMON /zet/     zets,deltaz,idz,ndeltaz,nidz,wzet
      COMMON /wv/      wavex,wavexs,wavey,waveys
      COMMON /chb/     chb
      COMMON /chbdx/   chbd1,chbd2
      COMMON /chndd/     cn_drdr,  ci_drdr,  e_drdr,
     &                   cn_nodr,  ci_nodr,  e_nodr,
     &                   cn_nono,  ci_nono,  e_nono,
     &                 uvcn_mimi,uvci_mimi,uve_mimi,
     &                  wcn_mimi, wci_mimi, we_mimi,
     &                 ttcn_mimi,ttci_mimi,tte_mimi
      COMMON /mltm/    mlt
      COMMON /pbot/    pbxnb,pcxnb,pbynb,pcynb,pbxob,pcxob,pbyob,pcyob
      COMMON /ptop/    pbxnt,pcxnt,pbynt,pcynt,pbxot,pcxot,pbyot,pcyot
      COMMON /prbc/      prbc1,  prbcn,
     &                 uvprbc1,uvprbcn,
     &                  wprbc1, wprbcn,
     &                 ttprbc1,ttprbcn
      COMMON /veluvw/  u,v,w
      COMMON /press/   p,rp
      COMMON /velruvw/ ru,rv,rw
      COMMON /adam/    h1,h2,h3
c      COMMON /large1g/ g
      COMMON /transp/  auxt1,auxt2
      COMMON /thrm0/   tt,rtt,htt,u_yz_plane,v_yz_plane,w_yz_plane
     &                ,u_2d_yz_plane,v_2d_yz_plane,w_2d_yz_plane
      COMMON /files/   ch_file,tt_file,format_mode,lag_file
      COMMON /velul/   ul
      COMMON /velvl/   vl
      COMMON /velwl/   wl
      COMMON /veldxul/ ddxul
      COMMON /veldyul/ ddyul
      COMMON /veldzul/ ddzul
      COMMON /spv1/    svx,svy
      COMMON /spv2/    msvx,msvy
      COMMON /spv3/    epssvx,epssvy
      COMMON /dep/     deposit,erosion
      COMMON /bc_v/    ubct,ubcb,ubctax,ubcbax,
     &                 vbct,vbcb,vbctax,vbcbax,
     &                 uvmbct1,uvmbcb1,uvmbctn,uvmbcbn
      COMMON /bc_w/    wbct,wbcb,wmbct1,wmbcb1,wmbctn,wmbcbn
      COMMON /bc_tt/   ttbct,ttbcb,ttmbct1,ttmbcb1,ttmbctn,ttmbcbn
      COMMON /grav/    gy,gz
      COMMON /kplus/   kts 
      COMMON /avgtime/ NT
      COMMON /avgvar/  um,vm,wm,ttm
     &                ,u2m,v2m,w2m,tt2m
     &                ,u3m,v3m,w3m,tt3m
     &                ,u4m,v4m,w4m,tt4m
     &                ,uvm,uwm,uttm
     &                ,vwm,vttm
     &                ,wttm
      COMMON /avgvar2/  uxm,uym,uzm
     &                 ,ux2m,uy2m,uz2m
     &                 ,ux3m,uy3m,uz3m
     &                 ,ux4m,uy4m,uz4m
      COMMON /avgvar3/  vxm,vym,vzm
     &                 ,vx2m,vy2m,vz2m
     &                 ,vx3m,vy3m,vz3m
     &                 ,vx4m,vy4m,vz4m
      COMMON /avgvar4/ wxm,wym,wzm
     &                ,wx2m,wy2m,wz2m
     &                ,wx3m,wy3m,wz3m
     &                ,wx4m,wy4m,wz4m
      COMMON /avgvar5/ ttxm,ttym,ttzm
     &                ,ttx2m,tty2m,ttz2m
     &                ,ttx3m,tty3m,ttz3m
     &                ,ttx4m,tty4m,ttz4m
      COMMON /avgvar6/ uespxm,uespym !,udspxm,udspym
     &                ,vespxm,vespym !,vdspxm,vdspym
     &                ,wespxm,wespym !,wdspxm,wdspym
      COMMON /avgvar7/ pwm,pm
     &                ,u2wm,v2wm
     &                ,uuzm,uwxm
     &                ,vvzm,vwym
     &                ,wwzm

ccccccccccccccccccccccc 2nd scalar field ccccccccccccccccccccccc
      COMMON /rparmt/  pr2,ra2,Tmax2
      COMMON /spv2/    msvx2,msvy2
      COMMON /spv3/    epssvx2,epssvy2
      COMMON /velul/   ul2
      COMMON /velvl/   vl2
      COMMON /velwl/   wl2
      COMMON /spv1/    svx2,svy2
      COMMON /press/   rp2
      COMMON /velruvw/ ru2,rv2,rw2
      COMMON /veldxul/ ddxul2

      COMMON /abg/     tt2ab,tt2at,tt2bb,tt2bt,tt2gb,tt2gt
      COMMON /chndd/   tt2cn_mimi,tt2ci_mimi,tt2e_mimi
      COMMON /prbc/    tt2prbc1,tt2prbcn
      COMMON /thrm0/   tt2,rtt2,htt2
      COMMON /files/   tt2_file
      COMMON /bc_tt/   tt2bct,tt2bcb,tt2mbct1,tt2mbcb1,tt2mbctn,tt2mbcbn
      COMMON /forcing/ forcing_tt2
      COMMON /ibmtt/  dtt2_c_dum_n,tt2p_tilde_n
ccccccccccccccccccccccc 2nd scalar field ccccccccccccccccccccccc

ccccccccccccccccccccccc Immersed boundary cccccccccccccccccccccc
      COMMON /ibm_tt/  tt_p,tt2_p
      COMMON /ibm_tt/  dttp_tildedn,dtt2p_tildedn
      COMMON /surf/    T_wall, Y_wall, T_inf, Y_inf 
      COMMON /surf/    z_o_avg, AF, diff_coef
      COMMON /surf/    n_l, num_x_o
      COMMON /surf/    deltaq, Rpillar
      COMMON /surf/    x01, y01, x02, y02

ccccccccccccccccccccccc Immersed boundary cccccccccccccccccccccc

c
cccccccccccccccccccccccccccccccccccccccccc GA cccccccccccccccccc
c     common variables used by the ibm
      COMMON /zet2/    xets,yets
      COMMON /bdld/    zplus,rplus
      

      COMMON /bdld/    z2d
      COMMON /bdld/    rf
      COMMON /bdld/    exner_rhs
      COMMON /forcing/ forcing_x,forcing_y,forcing_z,forcing_tt
      COMMON /delf/    ddf
      COMMON /sph/     irkk
      COMMON /lagr/    x_o,y_o,z_o,dv_l
      COMMON /lagr/    z_o_nm1, rho_fro_nm1
      COMMON /domainlag/  p_iw,p_ie,p_js,p_jn,p_kb,p_kt
      COMMON /domainlag/  p_iw_n,p_ie_n,p_js_n,p_jn_n,p_kb_n,p_kt_n
      COMMON /moving/  u_p,v_p,w_p
      COMMON /ibmtt/   chbd1_ibm,drag,lift
      COMMON /pibm/    dpx,dpy,dpz
cccccccccccccccccccccccccccccccccccccccccc GA

