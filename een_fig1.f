c     EE for ne fermions on a nx 1D chain
      
      program main 
      implicit complex*16 (z)
      implicit real*8(a-h,o-y)
      integer n,nx,ny,ne,ifail,ia,ib
      parameter (n=3003,n0=128,nx=14,ne=8,no=16384,ia1=64*n)
      integer a(n,nx),b(nx),jc(2),jd(2),an(no)
      integer a1(n0,nx/2),n1(n0),ai(nx),aj(nx)
      double precision dt,r,ec,phi1,phi2,flux,vx,vy,cc,dc,pi,ee1,cc1

      double precision eps,d(n),dr(n0),ux(n)
      double precision h(nx)
      double precision vvr(n,n)
      double precision ur(n0,n0)
      double precision wk11(ia1)

c     initialization      
      
      do 5 i2=1,nx
        do 5 i1=1,n
         a(i1,i2)=0
 5    continue

      do i1=1,no
         an(i1)=0
      end do

c Real space basis states with nx sites and up to ne electrons for example |100101110010101>
      
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

c     interaction contribution fo each basis state
      
      do i1=1,n
         ux(i1)=0.d0
         do i2=1,nx-1
            ux(i1)=ux(i1)+a(i1,i2)*a(i1,i2+1)
         end do
      end do

c Real space basis states with nx/2 sites and up to ne electrons for example |1000101> 
      
      l=0
      do 12 i1=0,1
       do 12 i2=0,1
        do 12 i3=0,1
         do 12 i4=0,1
          do 12 i5=0,1
           do 12 i6=0,1
            do 12 i7=0,1
c             do 12 i8=0,1

            l=l+1
            a1(l,1) = i1
            a1(l,2) = i2
            a1(l,3) = i3
            a1(l,4) = i4
            a1(l,5) = i5
            a1(l,6) = i6
            a1(l,7) = i7
c            a1(l,8) = i8

            n1(l)=i1+i2+i3+i4+i5+i6+i7
c            n1(l)=i1+i2+i3+i4+i5+i6+i7+i8
            
 12   continue 
      
c     number of samples
      
      nsamp = 1


c      pi=3.1415927d0

c     hopping amplitude
      dv=1.d0
c     interaction strength
      ud=.01d0

c     loop on samples (for averaging over disorder)
      
      do 113 id=1,nsamp

c     result file
         
      open(30,file='ferl14n8u001sse.x',
     .     form='FORMATTED',position='append') 

c     initilization of Hamiltoian matrix

      do 4 i1=1,n
       do 4 i2=1,n
        vvr(i1,i2)=0.d0
 4    continue

c     Interaction contribution to the Hamiltonian (always on diagonal)

      do i1=1,n
         vvr(i1,i1)=ux(i1)*ud
      end do

c     Hopping contribution to the Hamiltonian
      
      do 20 i1=1,n
       do 20 i2=i1+1,n
        j=0
         do 30 i4=1,nx
          j=a(i1,i4)*a(i2,i4) + j
          if (a(i1,i4).ne.a(i2,i4)) then
           if (a(i1,i4).eq.1) then
            ix1=i4
           else
            ix2=i4
           end if
          end if
 30     continue
        if (j.eq.(ne-1)) then
         if ((ix1+1).eq.ix2) then
            vvr(i1,i2) = -1.d0
         end if
        end if
 20   continue
      

      do 57 i3=1,n
       do 57 i4=1,i3-1
        vvr(i3,i4) = vvr(i4,i3)
 57   continue

c     Diagonalization of Hamiltonian using Lapack subroutine
      
      call DSYEV ('V','U',n,vvr,n,d,wk11,ia1,ifail)

c     Calculation of EE between halfs of the sample

c     Calculation of the reduced density matrix \rho_A for the jj eigenstate
      
      do jj=1,n/2

c        Initializtion of the reduced density matrix \rho_A
         
         do ii=1,n0
            do kk=1,n0
               ur(ii,kk)=0.d0
            end do
         end do

c        Trace over B region (sites larger than nx/2)
         
         do i=1,n0
            do j=1,n0
               if (n1(i).ne.n1(j)) then
                  ur(i,j)=0.d0
               else
                  do ii=1,nx/2
                     ai(ii)=a1(i,ii)
                     aj(ii)=a1(j,ii)
                  end do
                  do k=1,n0
c                  write (30,*) i,j,k
                     do ii=1,nx/2
                      ai(nx/2+ii)=a1(k,ii)
                      aj(nx/2+ii)=a1(k,ii)
                     end do
                     if (n1(k)+n1(i).eq.ne) then
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
                        ur(i,j)=ur(i,j)+vvr(i1,jj)*vvr(j1,jj)
                     end if
                  end do
               end if
            end do
         end do

c     Diagonalization of \rho_A using Lapack subroutine

         call DSYEV ('N','U',n0,ur,n0,dr,wk11,ia1,ifail)

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

