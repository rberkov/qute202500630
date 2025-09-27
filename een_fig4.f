      program main 
      implicit complex*16 (z)
      implicit real*8(a-h,o-y)
      integer n,nx,ny,ne,ifail,ia,ib
      parameter (n=12870,n0=256,nx=16,no=65536,ia1=64*n)
      integer a(n,nx),b(nx),jc(2),jd(2),an(no)
      integer a1(n0,nx/2),n1(n0),ai(nx),aj(nx)
      double precision dt,r,ec,phi1,phi2,flux,vx,vy,cc,dc,pi,ee1,cc1

      double precision eps,d(n),dr(n0)
      double precision h(nx)
      complex*16 zvc(n,n),zw(ia1)
      complex*16 jr(nx,nx,nx,nx)
      complex*16 uc(n0,n0)
      double precision rw(ia1)

c     initialization  

      do 5 i2=1,nx
        do 5 i1=1,n
         a(i1,i2)=0
 5    continue

      do i1=1,no
         an(i1)=0
      end do

c Real space basis states with nx sites and ne electrons for example |10010111001010101>

      l=0
      do 10 i1=1,nx
       do 10 i2=i1+1,nx
        do 10 i3=i2+1,nx
         do 10 i4=i3+1,nx
          do 10 i5=i4+1,nx
           do 10 i6=i5+1,nx
            do 10 i7=i6+1,nx
             do 10 i8=i7+1,nx
c              do 10 i9=i8+1,nx
c               do 10 i10=i9+1,nx
            l=l+1    
            a(l,i1) = 1
            a(l,i2) = 1
            a(l,i3) = 1
            a(l,i4) = 1
            a(l,i5) = 1
            a(l,i6) = 1
            a(l,i7) = 1
            a(l,i8) = 1
c            a(l,i9) = 1
c            a(l,i10) = 1
            
 10   continue 

c     translation table between binary enumeration of a state to consecutive enumeration
      
      do i1=1,n
         nl=1
         do i2=1,nx
            nl=nl+2**(i2-1)*a(i1,i2)
         end do
         an(nl)=i1
      end do
      
c     Real space basis states with nx/2 sites and up to ne electrons for example |1000101>

      l=0
      do 12 i1=0,1
       do 12 i2=0,1
        do 12 i3=0,1
         do 12 i4=0,1
          do 12 i5=0,1
           do 12 i6=0,1
            do 12 i7=0,1
             do 12 i8=0,1
c              do 12 i9=0,1
c               do 12 i10=0,1
            l=l+1
            a1(l,1) = i1
            a1(l,2) = i2
            a1(l,3) = i3
            a1(l,4) = i4
            a1(l,5) = i5
            a1(l,6) = i6
            a1(l,7) = i7
            a1(l,8) = i8
            n1(l)=i1+i2+i3+i4+i5+i6+i7+i8
c            n1(l)=i1+i2+i3+i4+i5+i6+i7+i8+i9+i10
            
 12   continue 

c     number of samples
      
      nsamp = 10


c      pi=3.1415927d0
      dv=1.d0

c     On site disorder - for canonical SYK always 0 
      dh =0.d0

c      initial sample
      
      ih=1091

c     loop on samples (for averaging over disorder)      
      
      do 113 id=1,nsamp

c     result file         
         
      open(30,file='syel16n8sse_1091.x',
     .     form='FORMATTED',position='append') 

c     temp file contains what sample is calculated now      
      
      open(31,file='syel16n8sse_1091.tmp',form='FORMATTED'
     .     ,position='append')

      write (31,*) id+ih

      close (31)

c    seed for random number generator
      
      iseed = 1111001+(id+ih)*2

      call srand(iseed)
      
c     for CSYK no on-site disorder
      
      do 51 i=1,nx
         h(i) = dh*ran(iseed) - dh/2.d0
 51   continue


c     initialization of complex Hamiltonian 
      
      do 4 i1=1,n
       do 4 i2=1,n
        zvc(i1,i2)=dcmplx(0.d0,0.d0)
 4    continue

c     initialization of complex 4 particle interaction 

      do 6 i1=1,nx
       do 6 i2=1,nx
        do 6 i3=1,nx
          do 6 i4=1,nx
            jr(i1,i2,i3,i4)=dcmplx(0.d0,0.d0)
 6    continue

      do 7 i1=1,nx
       do 7 i2=1,i1-1
        do 7 i3=1,nx
          do 7 i4=1,i3-1
             theta=2.d0*pi*rand()
             jr(i1,i2,i3,i4)=
     .       dcmplx(dv*dcos(theta),dv*dsin(theta))
     .            /(2.d0*nx)**1.5
 7       continue
         
c     Not relevant for CSYK
         
      do i1=1,n
       h1=0.d0
       do ix=1,nx
        h1=(dfloat(a(i1,ix))-.5d0)*h(ix)+h1
       end do
       zvc(i1,i1)=dcmplx(h1,0.d0)
      end do

c     setting the complex 4 particle interaction contribution to the Hamiltonian       
      
      do 20 i1=1,n
       do 20 i2=i1+1,n
        j=0
        j1=0
        j2=0
        do i=1,nx
          if (a(i1,i).ne.a(i2,i)) then
           j=j+1
           j1=j1+a(i1,i)
           j2=j2+a(i2,i)
          end if
        end do
        if (j.eq.4.and.j1.eq.2.and.j2.eq.2) then
         l1=1
         l2=1
         do i=1,nx
          b(i)=a(i2,i)
          if (a(i1,i).ne.a(i2,i)) then
           if (a(i2,i).eq.1) then
              jd(l1)=i
              l1=l1+1
           else
              jc(l2)=i
              l2=l2+1
           end if
          end if
         end  do
         nn=0
         do k=1,jd(1)-1
          nn=nn+b(k)
         end do
         b(jd(1))=0
         do k=1,jd(2)-1
           nn=nn+b(k)
         end do
         b(jd(2))=0
         do k=1,jc(1)-1
          nn=nn+b(k)
         end do
         b(jc(1))=1
         do k=1,jc(2)-1
           nn=nn+b(k)
         end do
         b(jc(2))=1
         zvc(i1,i2)=jr(jc(2),jc(1),jd(2),jd(1))*(-1.d0)**nn
        end if

 20   continue

      do 57 i3=1,n
       do 57 i4=1,i3-1
        zvc(i3,i4) = dconjg(zvc(i4,i3))
 57   continue


c     Diagonalization of Hamiltonian using Lapack subroutine      
      
      call ZHEEV ('V','U',n,zvc,n,d,zw,ia1,rw,ifail)

c     Calculation of EE between halfs of the sample

c     Calculation of the reduced density matrix \rho_A for the jj eigenstate
      
      do jj=1,n/2

c        Initializtion of the reduced density matrix \rho_A         
         
         do ii=1,n0
            do kk=1,n0
               uc(ii,kk)=(0.d0,0.d0)
            end do
         end do

c        Trace over B region (sites larger than nx/2)
         
         do i=1,n0
            do j=1,n0
               if (n1(i).ne.n1(j)) then
                  uc(i,j)=(0.d0,0.d0)
               else
                  do ii=1,nx/2
                     ai(ii)=a1(i,ii)
                     aj(ii)=a1(j,ii)
                  end do
                  do k=1,n0
                     do ii=1,nx/2
                      ai(nx/2+ii)=a1(k,ii)
                      aj(nx/2+ii)=a1(k,ii)
                     end do
                     if (n1(k)+n1(i).eq.nx/2) then
                        nl=1
                        do i1=1,nx
                           nl=nl+2**(i1-1)*ai(i1)
                        end do
                        nl1=1
                        do i1=1,nx
                           nl1=nl1+2**(i1-1)*aj(i1)
                        end do
                        i1=an(nl)
                        j1=an(nl1)
                        uc(i,j)=uc(i,j)+conjg(zvc(i1,jj))*zvc(j1,jj)
                     end if
                  end do
               end if
            end do
         end do

c     Diagonalization of \rho_A using Lapack subroutine         
         
         call ZHEEV ('N','U',n0,uc,n0,dr,zw,ia1,rw,ifail)

c     EE calculation         
         
      dd=0.d0
      do i=1,n0
         if (dr(i).gt.1.d-12) then
          dd=dd-dr(i)*dlog(dr(i))
         end if
       end do

       write (30,*) d(jj),dd

      end do

      write (30,*) ' '
      close (30)
      
 113  continue

      end

