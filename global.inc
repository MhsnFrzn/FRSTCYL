c  Global parameter file.
c
c  Number of grid points in each direction
c
      INTEGER nx,ny,nz
      PARAMETER (nx=300, ny=720 , nz=201)
cccccccccccccccccccccccccccccccccc GA ccccccccccc
c Specify the # of spheres(num_p) and approx.# of 
c Lagrangian points (n_ll)
      INTEGER num_p,n_ll,num_y_o,nxyp,nzp
      PARAMETER (nxyp=567, nzp=nz-6)
      PARAMETER (num_y_o=ny, num_p=8, n_ll=nx*ny)
cccccccccccccccccccccccccccccccccc GA
c
c  Define number of phases to average
c
      INTEGER np,dnp !recall that if np=1 => dnp=1
      PARAMETER (np=1,dnp=1)
c
c  Number of scalar fields
c 
      INTEGER m
      PARAMETER (m=1)
c
c  Derived parameters
c
      INTEGER nxl,nyl
      PARAMETER (nxl=3*nx/2, nyl=3*ny/2)
c
      INTEGER nxh,nyh,nxhp,nyhp,nxhm,cnxh,nzm,nzmm
      INTEGER nxlh,nylh,nxlhp,nylhp
      PARAMETER (nxh = nx/2, nyh = ny/2)
      PARAMETER (nxhp = nxh+1, nyhp = nyh+1)
      PARAMETER (nxhm = nxh-1, cnxh = nxh+2)
      PARAMETER (nxlh = nxl/2, nylh = nyl/2)
      PARAMETER (nxlhp = nxlh+1, nylhp = nylh+1)
      PARAMETER (nzm = nz-1, nzmm=nz-2)
c
      INTEGER nxlch,nxlchm,nxlchp,cnxlch,cnxlchm
      INTEGER nylch,nylchp
c
      PARAMETER (nxlch = nxh)
      PARAMETER (nxlchm = nxlch-1, nxlchp = nxlch+1)
      PARAMETER (cnxlch = nxl+2-nxlch, cnxlchm = cnxlch-1)
      PARAMETER (nylch = nyh, nylchp = nylch+1)
c
c  Non-final array subscript sizes (should avoid powers of two)
c
      INTEGER nx0,ny0,nz0,nzmm0,nxhp0,nyhp0,tnxhp0,tnyhp0
      INTEGER nxl0,nyl0,nxlhp0,nylhp0,tnxlhp0,tnylhp0
      PARAMETER (nx0=nx+1, ny0=ny+1, nz0=nz+1)
      PARAMETER (nzmm0=nzmm+1, nxhp0=nxhp, nyhp0=nyhp)
      PARAMETER (tnxhp0=2*nxhp0,tnyhp0=2*nyhp0)
      PARAMETER (nxl0=nxl+1, nyl0=nyl+1)
      PARAMETER (nxlhp0=nxlhp, nylhp0=nylhp)
      PARAMETER (tnxlhp0=2*nxlhp0,tnylhp0=2*nylhp0)
c
      INTEGER nxny,nx0y,nx0ylch,nx0y0,nx0y0z,nz0z
      INTEGER nyhp0x0,nylhp0xl0,nxl0yl
      PARAMETER (nxny=nx*ny,nx0y=nx0*ny,nx0ylch=nx0*nylch)
      PARAMETER (nx0y0=nx0*ny0,nx0y0z=nx0*ny0*nz,nz0z=nz0*nz)
      PARAMETER (nyhp0x0=nyhp0*nx0,nylhp0xl0=nylhp0*nxl0)
      PARAMETER (nxl0yl=nxl0*nyl)

