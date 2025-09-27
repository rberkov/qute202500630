c     Imbrie model for a chain of nx spins
      
      integer n,nx,ny,ne,ifail,ia,ib
      parameter (n=16384,n0=128,nx=14,nx0=7,ia1=64*n)
      integer a(n,nx), b(n0,nx0)
      double precision vvr(n,n),vj(nx)
      double precision ur(n0,n0)
      double precision dt,r,ec,phi1,phi2,flux,vx,vy,cc,dc,pi,ee1,cc1

      double precision eps,d(n),dr(n0)
      double precision wk11(64*n)
      double precision h(nx)

      
c     number of samples and starting sample
      
      nsamp = 100
      nsamp1 = 901

c      pi=3.1415927d0

c     disorder strength
      dh =3.d0

c     loop over samples
      
      do 113 id=nsamp1,nsamp+nsamp1-1

c     result file         
         
      open(30,file='imbl14pkheej04w3sse_0901.x',form='FORMATTED'
     .        ,position='append')

c     temp file - which sample is currently calculated
      
      open(31,file='imbl14pkheej04w3sse_0901.tmp',
     .     form='FORMATTED')
      write (31,*) id
      close (31)

c     seed for random number generation 
      
      iseed = 111100001+id*1002

      call srand(iseed)

c      initializing random magnetic fields h(i) 
c      and \gamma spin-spin interactions
      
      do 51 i=1,nx
        h(i) = dh*rand() - dh/2.d0
        vj(i)=.8d0+0.4*rand()
 51   continue

c     initilization of Hamiltoian matrix      

      do 4 i1=1,n
       do 4 i2=1,n
        vvr(i1,i2)=0.d0
 4    continue


      
      do 5 i2=1,nx
        do 5 i1=1,n
         a(i1,i2)=0
 5    continue

c  a represents real space basis states with nx spin 1/2 for example |100101110010101> 
c  b represents real space basis states with nx/2 spins 1/2 - for example |1000101>       
      
      l=0
      l0=0
      do 10 i1=0,1
       do 10 i2=0,1
        do 10 i3=0,1
         do 10 i4=0,1
          do 10 i5=0,1
           do 10 i6=0,1
            do 10 i7=0,1
            l0=l0+1
            b(l0,1)=i1
            b(l0,2)=i2
            b(l0,3)=i3
            b(l0,4)=i4
            b(l0,5)=i5
            b(l0,6)=i6
            b(l0,7)=i7
            do 10 i8=0,1
             do 10 i9=0,1
              do 10 i10=0,1
               do 10 i11=0,1
                do 10 i12=0,1
                 do 10 i13=0,1
                  do 10 i14=0,1
            l=l+1
            a(l,1) = i1
            a(l,2) = i2
            a(l,3) = i3
            a(l,4) = i4
            a(l,5) = i5
            a(l,6) = i6
            a(l,7) = i7
            a(l,8) = i8
            a(l,9) = i9
            a(l,10) = i10
            a(l,11) = i11
            a(l,12) = i12
            a(l,13) = i13
            a(l,14) = i14

c     magnetif field and \Gamma Interaction contribution to the Hamiltonian (always on diagonal)
            
            do jx=1,nx
             vvr(l,l)=(dfloat(a(l,jx))-.5d0)*h(jx)+vvr(l,l)
            end do
            do jx=1,nx-1
               vvr(l,l)=(dfloat(a(l,jx))-.5d0)*
     .                  (dfloat(a(l,jx+1))-.5d0)*vj(jx)+vvr(l,l)
            end do

 10   continue   

c     \sigma_x term contribution to the Hamiltonian
      
      do 20 i1=1,n
       do 20 i2=i1+1,n
        j=0
        do i=1,nx
          if(a(i1,i).eq.a(i2,i)) then
           j=j+1
          else
           jjx=j 
          end if
        end do
        if (nx-1.eq.j) then
            vvr(i1,i2)=.5d0
        end if
 20   continue

c     Diagonalization of Hamiltonian using Lapack subroutine
      
      call DSYEV ('V','U',n,vvr,n,d,wk11,ia1,ifail)

c     Calculation of EE between halfs of the sample

c     Calculation of the reduced density matrix \rho_A for the ii eigenstate 

      do ii=1,n

c     Trace over B region (spins # larger than nx/2)
         
      do i=1,n0
         i1=0
         do ix=1,nx0
            i1=i1+(2**(ix-1))*b(i,ix)
         end do
         do j=1,n0
           j1=0
           do ix=1,nx0
              j1=j1+(2**(ix-1))*b(j,ix)
           end do
           ur(i,j)=0.d0
           do k=1,n0
              k1=0
              do ix=1,nx0
                 k1=k1+(2**(ix+nx0-1))*b(k,ix)
              end do
              ii1=i1+k1+1
              jj1=j1+k1+1
              ur(i,j)=ur(i,j)+vvr(ii1,ii)*vvr(jj1,ii)
           end do
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

       write (30,*) d(ii),dd

c      print *,dd

      end do
      
      write (30,*) ' '
      close(30)      
 113  continue

      end

