      SUBROUTINE del_fn(l1,n1)
#include "common.inc"
c
c Dirac Delta Function based on Roma et al. (1999) and Hyungoo Lee
c
      INTEGER i,i1,i2,j,j1,j2,k,k1,k2,zref,l1,n1,i_t,j_t
      REAL*8 dum1(nx),dum2(ny),dum3(nz),dum333(nz),R,dzz,zminimum
      REAL*8 zarray(3),zarray1(2),hzp,hzm,hz
      REAL*8 ddf_control,dumc1,dumc2,dumc3
c
c Setting the limits of the Lagragian points domain to reduce
c computation
      i1=p_iw(l1,n1)
      i2=p_ie(l1,n1)
      j1=p_js(l1,n1)
      j2=p_jn(l1,n1)
      k1=p_kb(l1,n1)
      k2=p_kt(l1,n1)
c
c Calculating the dirac delta function in the x and y directions
      do i=i1,i2
         if (i .lt. 1) then
             R=((x_o(l1,n1)-xets(1))
     &              +dabs(xets(nx0)-xets(nx+i)))/deltax
             call FI_3point(R)
             dum1(nx+i)=FI/deltax
         elseif(i .gt. nx) then
             R=(dabs(xets(nx0)-x_o(l1,n1))
     &              +dabs(xets(i-nx)-xets(1)))/deltax
             call FI_3point(R)
             dum1(i-nx)=FI/deltax
         else
             R=dabs((xets(i)-x_o(l1,n1)))/deltax
             call FI_3point(R)
             dum1(i)=FI/deltax
         endif
      enddo
c
      do j=j1,j2
         if (j .lt. 1) then
             R=((y_o(l1,n1)-yets(1))
     &              +dabs(yets(ny0)-yets(ny+j)))/deltay
             call FI_3point(R)
             dum2(ny+j)=FI/deltay
         elseif(j .gt. ny) then
             R=(dabs(yets(ny0)-y_o(l1,n1))
     &              +dabs(yets(j-ny)-yets(1)))/deltay
             call FI_3point(R)
             dum2(j-ny)=FI/deltay
         else
             R=dabs((yets(j)-y_o(l1,n1)))/deltay
             call FI_3point(R)
             dum2(j)=FI/deltay
         endif
      enddo
c
ccc Calculating the local dialation factor (hz) in the z direction 
ccc ref. A. Pinelli_2010
c
c Finding the closest z cartesian point to the Lagrangian point
      zminimum=rlenz
      do k=1,nz
         dzz=dabs(zets(k)-z_o(l1,n1))
         if (dzz .lt. zminimum) then
            zref=k
            zminimum=dzz
         endif
      enddo
c
c defining the dialation factor hz
      if (zref .eq. 1) then
             hz=dabs(zets(zref+1)-zets(zref))/2.d0
             hzp=0
             hzm=0
      elseif (zref .eq. nz) then
             hz=dabs(zets(zref)-zets(zref-1))/2.d0
             hzp=0
             hzm=0
      elseif (zref .eq. 2) then
             zarray(1)=dabs(zets(2)-zets(1))
             zarray(2)=dabs(zets(3)-zets(2))
             zarray(3)=dabs(zets(4)-zets(3))
             hzp=max(zarray(1),zarray(2),zarray(3))
             hzm=min(zarray(1),zarray(2),zarray(3))
             hz=5.d0/6.d0*hzp+1.d0/6.d0*hzm
      else
             zarray(1)=dabs(zets(zref-1)-zets(zref-2))
             zarray(2)=dabs(zets(zref)-zets(zref-1))
             zarray(3)=dabs(zets(zref+1)-zets(zref))
             hzp=max(zarray(1),zarray(2),zarray(3))
             hzm=min(zarray(1),zarray(2),zarray(3))
             hz=5.d0/6.d0*hzp+1.d0/6.d0*hzm
      endif
c
c Calculating the dirac delta function in the z direction
      if (zref .eq. 1 .or. zref .eq. nz) then
         do k=k1,k2
            call FI_2point(z_o(l1,n1),k,zref)
            dum3(k)=FI
         enddo
      else
         do k=k1,k2
            R=dabs(zets(k)-z_o(l1,n1))/hz/1.0001d0
            call FI_3point(R)
            dum3(k)=FI
         enddo
      endif
c
c Modification of dum3 for non-uniform mesh
      do k=k1,k2
         dum33(k)=dum3(k)/hz/1.0001d0
c      write(117,*) l1,zref,dum3(k)
      enddo
c      write(172,*) zref
      if (zref .eq. 1 .or. zref .eq. nz) then
         do k=k1,k2
            dum333(k)=dum3(k)
         enddo
      else
c--------------------------------------------------
cc Correction for nonuniform z direction
         call modified_windows_fn(z_o(l1,n1),k1,k2,hz)
c--------------------------------------------------
         do k=k1,k2
            dum333(k)=(ibmb(1)+(zets(k)-z_o(l1,n1))*ibmb(2)
     &                +(zets(k)-z_o(l1,n1))**2.d0*ibmb(3))*dum33(k)
         enddo
      endif
c Combining the 1-D dirac delta functions to get the 3-D ddf
cC$OMP  PARALLEL DO DEFAULT(SHARED), PRIVATE(i1,i2,j1,j2,k1,k2
cC$OMP& ,i,j,k,i_t,j_t)
cC$OMP  PARALLEL DO DEFAULT(SHARED), PRIVATE(i,j,k,i_t,j_t)

      do i=i1,i2
         if(i .lt. 1) then
            i_t=nx+i
         elseif(i .gt. nx) then
            i_t=i-nx
         else
            i_t=i
         endif
         do j=j1,j2
            if(j .lt. 1) then
               j_t=ny+j
            elseif(j .gt. ny) then
               j_t=j-ny
            else
               j_t=j
            endif
            do k=k1,k2
               ddf_dum(i_t,j_t,k)=dum1(i_t)*dum2(j_t)*dum333(k)
            enddo
         enddo
      enddo
cC$OMP END PARALLEL DO
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
      return
      END SUBROUTINE del_fn
c
c
c
      SUBROUTINE FI_3point(R)
#include "common.inc"
c 3 point delta function by Peskin (2002)
      REAL*8 R
c
      if (R .le. 0.5d0) then
         FI=1.d0/3.d0*(1.d0+sqrt(-3.d0*R**2.d0+1.d0))
      elseif (R .le. 1.5d0) then
         FI=1.d0/6.0d0*(5.d0-3.d0*R-sqrt(-3.d0*(1.d0-R)**2.d0+1.d0))
      else
         FI=0.d0
      endif
      return
      END SUBROUTINE FI_3point
c
c
c
      SUBROUTINE FI_2point(z_0,k,reference)
#include "common.inc"
      REAL*8 z_0,sum
      INTEGER k,reference
c
      if (reference .eq. 1) then
         sum=dabs(z_0-zets(1))+dabs(z_0-zets(2))
         if (k .eq. 1) then
            FI=(zets(2)-z_0)/sum/(zets(2)-zets(1))*2.d0
         elseif (k .eq. 2) then
            FI=(z_0-zets(1))/sum/(zets(3)-zets(1))*2.d0
         else
            FI=0.d0
         endif
      endif
      if (reference .eq. nz) then
         sum=dabs(z_0-zets(nz))+dabs(z_0-zets(nzm)) 
         if (k .eq. nz) then
            FI=(zets(nzm)-z_0)/sum/(zets(nzm)-zets(nz))*2.d0
         elseif (k .eq. nzm) then
            FI=(z_0-zets(nz))/sum/(zets(nzmm)-zets(nz))*2.d0
         else
            FI=0.d0
         endif
      endif
c
      return
      END SUBROUTINE FI_2point
