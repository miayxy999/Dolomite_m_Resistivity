c  **************************  elecfem2d.f  ********************************
c  BACKGROUND

c  This program solves Laplace's equation in a random conducting
c  material using the finite element method.  Each pixel in the 2-D digital
c  image is a square bi-linear finite element,  having its own conductivity
c  tensor.  Periodic boundary conditions are maintained.
c  In the comments below, (USER) means that this is a section of code
c  that the user might have to change for his particular problem.
c  Therefore the user is encouraged to search for this string.

c  PROBLEM AND VARIABLE DEFINITION

c  The problem being solved is the minimization of the energy
c  1/2 uAu + b u + C, where A is the Hessian matrix composed of the 
c  stiffness matrices (dk) for each pixel/element, b is a constant vector
c  and C is a constant that are determined by the applied field and 
c  the periodic boundary conditions, and u is a vector of all the voltages.
c  The method used is the conjugate gradient relaxation algorithm.
c  Other variables are:  gb is the gradient = Au+b, h and Ah are 
c  auxiliary variables used in the conjugate gradient algorithm (in dembx),
c  dk(n,i,j) is the stiffness matrix of the n'th phase, sigma(n,i,j) is
c  the conductivity tensor of the n'th phase, pix is a vector that gives
c  the phase label of each pixel, ib is a matrix that gives the labels of
c  the 9 (counting itself) neighbors of a given node, prob is the area 
c  fractions of the various phases, and currx, curry, currz are the 
c  area averaged total currents in the x, y, and z directions.

c  DIMENSIONS

c  The vectors u,gb,b,h, and Ah are dimensioned to be the system size,
c  ns=nx*ny, where the digital image of the microstructure considered
c  is a rectangle ( nx x ny) in size. The arrays ib and pix are also
c  dimensioned to the system size. The array ib has 9 components, 
c  for the 9 neighbors of a node.
c  Note that the program is set up at present to have at most 100 
c  different phases.  This can easily be changed, simply by changing 
c  the dimensions of dk, prob, and sigma. Nphase gives the number of
c  phases being considered.
c  All arrays are passed between subroutines using simple common statements.

c  STRONGLY SUGGESTED:  READ THE MANUAL BEFORE USING PROGRAM!!

c  (USER)  Change these dimensions and in other subroutines at the
c  same time.  For example, search and replace all occurrences throughout
c  the program of "(400" by "(1600", to go from a 20 x 20 system to
c  a 40 x 40 system.
	double precision u(3125000),gb(3125000),b(3125000),dk(100,4,4)
      double precision h(3125000),Ah(3125000)
	double precision sigma(100,2,2), prob(100)
      double precision currx, curry, ex, ey
      double precision utot
	integer in(9),jn(9),kn(9)
	integer*4 ib(3125000,9)
      integer*2 pix(3125000)


	common/list1/currx,curry
	common/list2/ex,ey
	common/list3/ib
      common/list4/pix
	common/list5/dk,b,C
	common/list6/u
	common/list7/gb
	common/list8/sigma
	common/list9/h,Ah

c  (USER)  Unit 9 is the microstructure input file, unit 7 is
c  the results output file.
        open(9,file='microstructure323_jizhi35s1.dat')
        open(7,file='outputfile323_jizhi35s1.out')

c  (USER) nx,ny give the size of the lattice
      nx=2500
      ny=1250
c ns=total number of sites
        ns=nx*ny
        write(7,9010) nx,ny,ns
9010    format('nx= ',i8,' ny= ',i8,' ns = ',i8)

c  (USER) nphase is the number of phases being considered in the problem.
c  The values of pix(m) will run from 1 to nphase.
      nphase=3

c  (USER) gtest is the stopping criterion, compared to gg=gb*gb.  
c  gtest=abc*ns, so that when gg < gtest, that average value per pixel
c  of gb is less than sqrt(abc).gtest=1.e-16*ns
      gtest=1.d-16*ns

c  Construct the neighbor table, ib(m,n)

c  First construct 9 neighbor table in terms of delta i, delta j
c  (See Table 3 in manual)
      in(1)=0
      in(2)=1
      in(3)=1
      in(4)=1
      in(5)=0
      in(6)=-1
      in(7)=-1
      in(8)=-1
      in(9)=0

      jn(1)=1
      jn(2)=1
      jn(3)=0
      jn(4)=-1
      jn(5)=-1
      jn(6)=-1
      jn(7)=0
      jn(8)=1
      jn(9)=0

c  Now construct neighbor table according to 1-d labels
c  Matrix ib(m,n) gives the 1-d label of the n'th neighbor (n=1,9) of
c  the node labelled m.
      do 1020 j=1,ny
      do 1021 i=1,nx
      m=nx*(j-1)+i
      do 1004 n=1,9
      i1=i+in(n)
      j1=j+jn(n)
      if(i1.lt.1) i1=i1+nx
      if(i1.gt.nx) i1=i1-nx
      if(j1.lt.1) j1=j1+ny
      if(j1.gt.ny) j1=j1-ny
      m1=nx*(j1-1)+i1
      ib(m,n)=m1
1004  continue
1021  continue 
1020  continue 

c  Compute the electrical conductivity of each microstructure
c  (USER) npoints is the number of microstructures to use
        npoints=1
        do 8000 micro=1,npoints

c  Read in a microstructure in subroutine ppixel, and set up pix(m) 
c  with the appropriate phase assignments.
        call ppixel(nx,ny,ns,nphase)
c Count and output the area fractions of the different phases
        call assig(ns,nphase,prob)
	do 805 i=1,nphase
	write(7,*) 'Area fraction of phase ',i,' = ',prob(i)
805     continue


c  (USER) sigma(100,2,2) is the electrical conductivity tensor of each phase
c  The user can make the value of sigma to be different for each 
c  phase of the microstructure if so desired, up to the limit of the
c  dimension of sigma.
      sigma(1,1,1)=18.717d0
      sigma(1,2,2)=18.717d0
      sigma(1,1,2)=0.0d0
      sigma(1,2,1)=sigma(1,1,2)

      sigma(2,1,1)=4.7d-4
      sigma(2,2,2)=4.7d-4
      sigma(2,1,2)=0.0d0
      sigma(2,2,1)=sigma(2,1,2)

      sigma(3,1,1)=1.0d-6
      sigma(3,2,2)=1.0d-6
      sigma(3,1,2)=0.0d0
      sigma(3,2,1)=sigma(3,1,2)

c  write out the phase electrical conductivity tensors
      do 11 i=1,nphase 
      write(7,*) 'Phase ',i,' conductivity tensor is:'
      write(7,*) sigma(i,1,1),sigma(i,1,2)
      write(7,*) sigma(i,2,1),sigma(i,2,2)
11    continue

c  (USER) Set the applied electric field.
          ex=1.0
          ey=1.0
      write(7,*) 'Applied field components:'
      write(7,*) 'ex = ',ex,' ey = ',ey

c  Set up the finite element "stiffness" matrices and the Constant and
c  vector required for the energy

	call femat(nx,ny,ns,nphase)

c  Apply a homogeneous macroscopic electric field as the initial condition
	do 1050 j=1,ny
        do 1051 i=1,nx
		m=nx*(j-1)+i
		x=float(i-1)
		y=float(j-1)
		u(m)=-x*ex-y*ey
1051	continue          
1050	continue

c  Relaxation Loop
c  (USER) kmax is the maximum number of times dembx will be called, with
c  ldemb conjugate gradient steps done during each call. The total
c  number of conjugate gradient cycles allowed for a given conductivity
c  computation is kmax*ldemb.
        kmax=4000
        ldemb=100
        ltot=0

c  Call energy to get initial energy and initial gradient
        call energy(nx,ny,ns,utot)
c  gg is the norm squared of the gradient (gg=gb*gb)
 	gg=0.0
        do 100 m=1,ns
        gg=gg+gb(m)*gb(m)
100     continue
	write(7,*) 'Initial energy = ',utot,'gg = ',gg
        call flush(7)

        do 5000 kkk=1,kmax
c  Call dembx to go into conjugate gradient solver
        call dembx(ns,Lstep,gg,dk,gtest,ldemb,kkk)
        ltot=ltot+Lstep
c  Call energy to compute energy after dembx call. If gg < gtest, this
c  will be the final energy. If gg is still larger than gtest, then this
c  will give an intermediate energy with which to check how the relaxation 
c  process is coming along.
        call energy(nx,ny,ns,utot)
	write(7,*) 'Energy = ',utot,'gg = ',gg
	write(7,*) ltot, ' conj. grad. steps'
        if(gg.lt.gtest) goto 444

c  Compute and output the currents as an additional aid to judge how
c  the relaxation process is progressing.
       call current(nx,ny,ns)
c Output currents
        write(7,*)
	write(7, *) ' Current in x direction = ',currx
	write(7, *) ' Current in y direction = ',curry
        call flush(7)

5000   continue

444    call current(nx,ny,ns)

c Output currents
        write(7,*)
	write(7, *) ' Current in x direction = ',currx
	write(7, *) ' Current in y direction = ',curry
        call flush(7)

8000    continue

        end

c  Subroutine that sets up the stiffness matrices, linear term in 
c  voltages, and constant term C that appear in the total energy due 
c  to the periodic boundary conditions.

      subroutine femat(nx,ny,ns,nphase)
      double precision dk(100,4,4),xn(8),b(3125000),C
      double precision dndx(4),dndy(4),dndz(4)
      double precision g(3,3),sigma(100,2,2)
      double precision es(2,4),sum,x,y
      double precision ex, ey
      double precision utot
      integer is(4)
      integer*4 ib(3125000,9)
      integer*2 pix(3125000)

	common/list2/ex,ey
	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C
	common/list8/sigma

c  initialize stiffness matrices
      do 40 m=1,nphase
      do 41 j=1,4
      do 42 i=1,4
      dk(m,i,j)=0.0
42    continue
41    continue
40    continue

c  set up Simpson's rule integration weight vector
      do 30 j=1,3
      do 31 i=1,3
      nm=0
      if(i.eq.2) nm=nm+1
      if(j.eq.2) nm=nm+1
      g(i,j)=4.0**nm
31    continue
30    continue

c  loop over the nphase kinds of pixels and Simpson's rule quadrature
c  points in order to compute the stiffness matrices.  Stiffness matrices
c  of bilinear finite elements are quadratic in x and y, so that
c  Simpson's rule quadrature is exact.
      do 4000 ijk=1,nphase
      do 3001 j=1,3
      do 3000 i=1,3
      x=float(i-1)/2.0
      y=float(j-1)/2.0
c  dndx means the negative derivative with respect to x of the shape 
c  matrix N (see manual, Sec. 2.2), dndy is similar.
      dndx(1)=-(1.0-y)
      dndx(2)=(1.0-y)
      dndx(3)=y
      dndx(4)=-y
      dndy(1)=-(1.0-x)
      dndy(2)=-x
      dndy(3)=x
      dndy(4)=(1.0-x)
c  now build electric field matrix
      do 2799 n1=1,2
      do 2799 n2=1,4
      es(n1,n2)=0.0
2799  continue
      do 2797 n=1,4
      es(1,n)=dndx(n)
      es(2,n)=dndy(n)
2797  continue
c  now do matrix multiply to determine value at (x,y), multiply by
c  proper weight, and sum into dk, the stiffness matrix
      do 900 ii=1,4
      do 900 jj=1,4
c  Define sum over field matrices and conductivity tensor that defines
c  the stiffness matrix.
      sum=0.0d0
      do 890 kk=1,2
      do 890 ll=1,2
      sum=sum+es(kk,ii)*sigma(ijk,kk,ll)*es(ll,jj)
890   continue
      dk(ijk,ii,jj)=dk(ijk,ii,jj)+g(i,j)*sum/36.0d0
900   continue
3000  continue
3001  continue
4000  continue

c  Set up vector for linear term, b, and constant term, C, 
c  in the electrical energy.  This is done using the stiffness matrices,
c  and the periodic terms in the applied field that come in at the boundary
c  pixels via the periodic boundary conditions and the condition that
c  an applied macroscopic field exists (see Sec. 2.2 in manual).

      do 5000 m=1,ns
      b(m)=0.0
5000  continue

c  For all cases, correspondence between 1-8 finite element node labels
c  and 1-27 neighbor labels is:  1:ib(m,9),2:ib(m,3),3:ib(m,2),
c  4:ib(m,1) (see Table 4 in manual)
      is(1)=9
      is(2)=3
      is(3)=2
      is(4)=1

      C=0.0
c  x=nx face
      i=nx
      do 2001 i4=1,4
      xn(i4)=0.0
      if(i4.eq.2.or.i4.eq.3) then
      xn(i4)=-ex*nx
      end if
2001  continue
      do 2000 j=1,ny-1
      m=j*nx
      do 1900 mm=1,4
      sum=0.0
      do 1899 m4=1,4
      sum=sum+xn(m4)*dk(pix(m),m4,mm)
      C=C+0.5*xn(m4)*dk(pix(m),m4,mm)*xn(mm)
1899  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1900  continue
2000  continue
c  y=ny face
      j=ny
      do 2011 i4=1,4
      xn(i4)=0.0
      if(i4.eq.3.or.i4.eq.4) then
      xn(i4)=-ey*ny
      end if
2011  continue
      do 2010 i=1,nx-1
      m=nx*(ny-1)+i
      do 1901 mm=1,4
      sum=0.0
      do 2099 m4=1,4
      sum=sum+xn(m4)*dk(pix(m),m4,mm)
      C=C+0.5*xn(m4)*dk(pix(m),m4,mm)*xn(mm)
2099  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1901  continue
2010  continue
c  x=nx y=ny corner 
      i=nx
      j=ny
      do 2061 i4=1,4
      xn(i4)=0.0
      if(i4.eq.2) then
      xn(i4)=-ex*nx
      end if
      if(i4.eq.4) then
      xn(i4)=-ey*ny
      end if
      if(i4.eq.3) then
      xn(i4)=-ex*nx-ey*ny
      end if
2061  continue
      m=nx*ny
      do 1906 mm=1,4
      sum=0.0
      do 2059 m4=1,4
      sum=sum+xn(m4)*dk(pix(m),m4,mm)
      C=C+0.5*xn(m4)*dk(pix(m),m4,mm)*xn(mm)
2059  continue
      b(ib(m,is(mm)))=b(ib(m,is(mm)))+sum
1906  continue

      return
      end

c  Subroutine computes the total energy, utot, and gradient, gb

      subroutine energy(nx,ny,ns,utot)
	double precision u(3125000),gb(3125000)
	double precision b(3125000),C
	double precision dk(100,4,4)
	double precision utot
      double precision ex, ey
 	integer*4 ib(3125000,9)
        integer*2 pix(3125000)
	
	common/list2/ex,ey
	common/list3/ib
        common/list4/pix
	common/list5/dk,b,C
	common/list6/u
	common/list7/gb

	do 2090 m=1,ns
	gb(m)=0.0
2090	continue

c  Energy loop. Do global matrix multiply via small stiffness matrices,
c  gb=Au + b.  The long statement below correctly brings in all the terms
c  from the global matrix A using only the small stiffness matrices.
      do 3000 m=1,ns
	gb(m)=gb(m)+u(ib(m,1))*(dk(pix(ib(m,9)),1,4)+
     &	dk(pix(ib(m,7)),2,3))+
     &  u(ib(m,2))*dk(pix(ib(m,9)),1,3)+
     $  u(ib(m,3))*
     &  (dk(pix(ib(m,9)),1,2)+
     &  dk(pix(ib(m,5)),4,3))+
     &  u(ib(m,4))*dk(pix(ib(m,5)),4,2)+
     &  u(ib(m,5))*(dk(pix(ib(m,5)),4,1)+dk(pix(ib(m,6)),3,2))+
     &  u(ib(m,6))*dk(pix(ib(m,6)),3,1)+
     &  u(ib(m,7))*
     &  (dk(pix(ib(m,6)),3,4)+dk(pix(ib(m,7)),2,1))+
     &  u(ib(m,8))*dk(pix(ib(m,7)),2,4)+
     &  u(ib(m,9))*(dk(pix(ib(m,9)),1,1)+
     &  dk(pix(ib(m,7)),2,2)+dk(pix(ib(m,6)),3,3)+
     &  dk(pix(ib(m,5)),4,4))
3000  continue

	utot=0.0
	do 3100 m=1,ns
	utot=utot+0.5*u(m)*gb(m)+b(m)*u(m)
	gb(m)=gb(m)+b(m)
3100	continue

	utot=utot+C

        return
        end           	       

c    Subroutine that carries out the conjugate gradient relaxation process

      subroutine dembx(ns,Lstep,gg,dk,gtest,ldemb,kkk)
      double precision u(3125000),gb(3125000),dk(100,4,4)    
      double precision Ah(3125000),h(3125000),B,lambda,gamma
      double precision utot
      integer*4 ib(3125000,9)
      integer*2 pix(3125000)

	common/list3/ib
        common/list4/pix
	common/list6/u
	common/list7/gb
	common/list9/h,Ah
      
c  Initialize the conjugate direction vector on first call to dembx only.
c  For calls to dembx after the first, we want to continue using the value 
c  of h determined in the previous call.  Of course, if npoints is greater
c  than 1, then this initialization step will be run each time a new
c  microstructure is used, as kkk will be reset to 1 every time the counter
c  micro is increased.
      if(kkk.eq.1) then
      do 50 m=1,ns
      h(m)=gb(m) 
50    continue                                                          
      end if
c  Lstep counts the number of conjugate gradient steps taken in each call
c  to dembx.  
      Lstep=0                                                              

c     Conjugate gradient loop

      do 800 ijk=1,ldemb
      Lstep=Lstep+1
          
      do 290 m=1,ns
      Ah(m)=0.0
290   continue
c  Do global matrix multiply via small stiffness matrices, Ah = A * h.
c  The long statement below correctly brings in all the terms from
c  the global matrix A using only the small stiffness matrices.
      do 400 m=1,ns
	Ah(m)=Ah(m)+h(ib(m,1))*(dk(pix(ib(m,9)),1,4)+
     &	dk(pix(ib(m,7)),2,3))+
     &  h(ib(m,2))*dk(pix(ib(m,9)),1,3)+
     &  h(ib(m,3))*
     &  (dk(pix(ib(m,9)),1,2)+
     &  dk(pix(ib(m,5)),4,3))+
     &  h(ib(m,4))*dk(pix(ib(m,5)),4,2)+
     &  h(ib(m,5))*(dk(pix(ib(m,5)),4,1)+dk(pix(ib(m,6)),3,2))+
     &  h(ib(m,6))*dk(pix(ib(m,6)),3,1)+
     &  h(ib(m,7))*
     &  (dk(pix(ib(m,6)),3,4)+dk(pix(ib(m,7)),2,1))+
     &  h(ib(m,8))*dk(pix(ib(m,7)),2,4)+
     &  h(ib(m,9))*(dk(pix(ib(m,9)),1,1)+
     &  dk(pix(ib(m,7)),2,2)+dk(pix(ib(m,6)),3,3)+
     &  dk(pix(ib(m,5)),4,4))
400   continue

      hAh=0.0                                                            
      do 530 m=1,ns                                                    
      hAh=hAh+h(m)*Ah(m)                                                 
530   continue                                                          

      lambda=gg/hAh
      do 540 m=1,ns                                                    
      u(m)=u(m)-lambda*h(m)                                        
      gb(m)=gb(m)-lambda*Ah(m)                                        
540   continue
                                                                
      gglast=gg
      gg=0.0
      do 550 m=1,ns                                                     
      gg=gg+gb(m)*gb(m)       
550   continue
      if(gg.le.gtest) goto 1000

      gamma=gg/gglast
      do 570 m=1,ns                                                   
      h(m)=gb(m)+gamma*h(m)                                       
570   continue                                                        

800   continue

1000  continue                                                        

      return                                                           
      end        

c Subroutine that computes average current in two directions

      subroutine current(nx,ny,ns)

      double precision af(2,4)
      double precision u(3125000),uu(4)
      double precision sigma(100,2,2)
      double precision currx, curry, ex, ey, cur1, cur2
      double precision utot
      integer*4 ib(3125000,9)
      integer*2 pix(3125000)

	common/list1/currx,curry
	common/list2/ex,ey
	common/list3/ib
        common/list4/pix
	common/list6/u
	common/list8/sigma

c  af is the average field matrix, average field in a pixel is af*u(pixel).
c  af relates the nodal voltages to the average field in the pixel.

c Set up single element average field matrix

      af(1,1)=0.5
      af(1,2)=-0.5
      af(1,3)=-0.5
      af(1,4)=0.5
      af(2,1)=0.5
      af(2,2)=0.5
      af(2,3)=-0.5
      af(2,4)=-0.5

c  now compute current in each pixel
      currx=0.0
      curry=0.0
c  compute average field in each pixel
      do 470 j=1,ny
      do 471 i=1,nx
      m=(j-1)*nx+i
c  load in elements of 4-vector using pd. bd. conds.
      uu(1)=u(m) 
      uu(2)=u(ib(m,3))
      uu(3)=u(ib(m,2))
      uu(4)=u(ib(m,1))
c  Correct for periodic boundary conditions, some voltages are wrong
c  for a pixel on a periodic boundary. Since they come from an opposite
c  face, need to put in applied fields to correct them.
      if(i.eq.nx) then 
      uu(2)=uu(2)-ex*nx
      uu(3)=uu(3)-ex*nx
      end if
      if(j.eq.ny) then
      uu(3)=uu(3)-ey*ny
      uu(4)=uu(4)-ey*ny
      end if
c  cur1 and cur2 are the local currents averaged over the pixel
      cur1=0.0
      cur2=0.0
      do 465 n=1,4
      do 465 nn=1,2
      cur1=cur1+sigma(pix(m),1,nn)*af(nn,n)*uu(n)
      cur2=cur2+sigma(pix(m),2,nn)*af(nn,n)*uu(n)
465   continue
c  sum into the global average currents
      currx=currx+cur1
      curry=curry+cur2
471   continue      
470   continue

c Area average currents
      currx=currx/dble(ns)
      curry=curry/dble(ns)

      return
      end

c  Subroutine that counts phase area fractions

      subroutine assig(ns,nphase,prob)

      integer ns,nphase
      integer*2 pix(3125000)
      double precision prob(100)
      integer*4 total1, total2, total3
      
      

      common/list4/pix
      total1 = 0
      total2 = 0
      total3 = 0
	do 90 i=1,nphase
	prob(i)=0.0
90      continue

      do 100 m=1,ns
        if(pix(m).eq.1) then
      total1 = total1 + 1
      endif
        if(pix(m).eq.2) then
      total2 = total2 + 1
      endif
        if(pix(m).eq.3) then
      total3 = total3 + 1
      endif
c	do 101 i=1,nphase
c        if(pix(m).eq.i) then
c	prob(i)=prob(i)+1
c	endif
101   continue      
100   continue

c	do 110 i=1,nphase
c	prob(i)=dble(prob(i))/dble(ns)
110     continue
      write(7,*) 'Total pixels of phase 1 = ', total1
      write(7,*) 'Total pixels of phase 2 = ', total2
      write(7,*) 'Total pixels of phase 3 = ', total3
      write(7,*) 'Sum = ', total1 + total2 + total3
      prob(1) = dble(total1) / dble(total1 + total2 + total3)
      prob(2) = dble(total2) / dble(total1 + total2 + total3)
      prob(3) = dble(total3) / dble(total1 + total2 + total3)
      return
      end

c  Subroutine that sets up microstructural image

      subroutine ppixel(nx,ny,ns,nphase)
      integer*2 pix(3125000)
      common/list4/pix

c  (USER)  If you want to set up a test image inside the program, instead
c  of reading it in from a file, this should be done inside this
c  subroutine.

      do 100 j=1,ny
      do 200 i=1,nx 
      m=nx*(j-1)+i
      read(9,*) pix(m)
200   continue   
100   continue

c  Check for wrong phase labels--less than 1 or greater than nphase
       do 500 m=1,ns
       if(pix(m).lt.1) then
        write(7,*) 'Phase label in pix < 1--error at ',m
       end if
       if(pix(m).gt.nphase) then
        write(7,*) 'Phase label in pix > nphase--error at ',m
       end if
500    continue

      return
      end
