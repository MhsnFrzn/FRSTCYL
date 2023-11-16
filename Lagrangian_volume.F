      SUBROUTINE lagrangian_volume(n1,ph,tht)
#include "common.inc"
c assign discrete volumes for the Lagrangian points
c Akiki and Balchandar, Vector Spherical Harmonics
c
      INTEGER   n1,neq,i,j,k
      REAL*8    ph(n_l(n1)),tht(n_l(n1))
      PARAMETER (neq=9)
      REAL*8    AA(neq,n_l(n1)),Source(neq),BB(n_l(n1),neq),
     &          AB(neq,neq),Areadum(neq),Area(n_l(n1))
      REAL*8    Pdum,PDdum,PP(4,3,n_l(n1))
      INTEGER   lwork,i22,i33,ipiv1(neq)
      PARAMETER (lwork=neq*(neq+4))
      REAL*8    work(lwork),thickness
c
      do i=1,neq
         do j=1,n_l(n1)
            AA(i,j)=0.d0
         enddo
         Source(i)=0.d0
         Areadum(i)=1.d0
      enddo
c
      do j=1,n_l(n1)
         AA(1,j)=dsin(tht(j))*dcos(ph(j))
         AA(2,j)=dsin(tht(j))*dsin(ph(j))
         AA(3,j)=dcos(tht(j))
         AA(4,j)=-dsin(tht(j))**2.d0*dcos(ph(j))**2.d0
         AA(5,j)=-dsin(tht(j))*dcos(ph(j))*dcos(tht(j))
         AA(6,j)=-dsin(tht(j))**2.d0*dsin(ph(j))**2.d0
         AA(7,j)=-dcos(tht(j))**2.d0*dcos(ph(j))**2.d0
     &           -(dsin(ph(j))**2.d0)
c         AA(8,j)=(-dcos(tht(j))**2.d0*dcos(ph(j))**2.d0
c     &           -dsin(ph(j))**2.d0)*rz_l(j,n1)
c     &    -(dcos(tht(j))*dcos(ph(j))*dsin(tht(j)))*rx_l(j,n1)
         AA(8,j)=-dcos(tht(j))**2.d0*dsin(ph(j))*dcos(ph(j))
     &           +dcos(ph(j))*dsin(ph(j))
         AA(9,j)=dcos(tht(j))*dsin(ph(j))*dsin(tht(j))
      enddo
c
      Source(4)=-4.d0/3.d0*pi*r_c**2.d0
      Source(6)=-4.d0/3.d0*pi*r_c**2.d0
      Source(7)=-8.d0/3.d0*pi*r_c**2.d0
c
      do i=1,3
         do j=1,i
            do k=1,n_l(n1)
               call plgndr(i-1,j-1,dcos(tht(k)),Pdum)
               PP(i,j,k)= Pdum
            enddo
         enddo
      enddo
c
      do i=1,n_l(n1)
         BB(i,1) = PP(1,1,i)
         BB(i,2) = PP(2,1,i)
         BB(i,3) = PP(2,2,i)*dcos(ph(i))
         BB(i,4) = PP(2,2,i)*dsin(ph(i))
         BB(i,5) = PP(3,1,i)
         BB(i,6) = PP(3,2,i)*dcos(ph(i))
         BB(i,7) = PP(3,2,i)*dsin(ph(i))
         BB(i,8) = PP(3,3,i)*dcos(2.d0*ph(i))
         BB(i,9) = PP(3,3,i)*dsin(2.d0*ph(i))
      enddo
c
      call dgemm('N','N',neq,neq,n_l(n1),1.d0,AA,neq,
     &                   BB,n_l(n1),0.d0,AB,neq)
      call dgetrf(neq,neq,AB,neq,ipiv1,i22)
      call dgetri(neq,AB,neq,ipiv1,work,lwork,i33)
      call dgemv('N',neq,neq,1.d0,AB,neq,Source,1,0.d0,Areadum,1)
      call dgemv('N',n_l(n1),neq,1.d0,BB,n_l(n1),Areadum,1,0.d0,Area,1)
c
      thickness=dsqrt(4.d0*pi*r_c**2.d0/n_l(n1))
      do i=1,n_l(n1)
         dv_l(i,n1)=Area(i)*thickness
      enddo
      do i=1,n_l(n1)
         write(183,*) dv_l(i,n1)
      enddo
c
      End
c
      SUBROUTINE plgndr(l,m,x,PPdum)
      INTEGER l,m
      REAL*8 x
      INTEGER i,ll
      REAL*8 fact,pll,pmm,pmmp1,somx2,PPdum
      if(m.lt.0.or.m.gt.l.or.abs(x).gt.1.)pause
c     *'bad arguments in plgndr'
      pmm=1.
      if(m.gt.0) then
        somx2=sqrt((1.-x)*(1.+x))
        fact=1.
        do 11 i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.
11      continue
      endif
      if(l.eq.m) then
c        plgndr=pmm
        PPdum=pmm
      else
        pmmp1=x*(2*m+1)*pmm
        if(l.eq.m+1) then
c          plgndr=pmmp1
          PPdum=pmmp1
        else
          do 12 ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
12        continue
c          plgndr=pll
          PPdum=pll
        endif
      endif
      return
      END SUBROUTINE plgndr
