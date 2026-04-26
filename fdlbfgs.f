c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine fdlbfgs 
     & (n,m,x,f,g,dfu,qpr,eps,lout,lerr,maxit,iter,warn)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c-----This routine is a simple driver for the lbfgs routine, a limited 
c     memory lbfgs method for large scale optimizations developed by 
c     Jorge Nocedal (1990). It solves the unconstrained minimization
c     problem  MIN F(x1,x2,...,xN) using the limited memory BFGS method.
c     The routine is especially effective on problems involving a large 
c     number of variables. In a typical iteration of this method an 
c     approximation Hk to the inverse of the Hessian is obtained by 
c     applying M BFGS updates to a diagonal matrix Hk0, using informa-
c     tion from the previous M steps. The user specifies the number M, 
c     which determines the amount of storage required by the routine. 
c     The user may also provide the diagonal matrices Hk0 if not satis-
c     fied with the default choice. The algorithm is described in 
c     "On the limited memory BFGS method for large scale optimization", 
c     by D. Liu and J. Nocedal, Mathematical Programming B 45 (1989) 
c     503-528.
c 
c     The user is required to calculate the function value F and its
c     gradient G. In order to allow the user complete control over
c     these computations, reverse  communication is used. The lbfgs
c     routine must be called repeatedly under the control of the 
c     parameter iflag. 
c
c     The steplength is determined at each iteration by means of the
c     line search routine MCVSRCH, which is a slight modification of
c     the routine CSRCH written by More and Thuente.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     DESCRIPTION OF THE INPUT PARAMETERS.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c 
c-----N       
c
c     INTEGER variable that must be set by the user to the number of 
c     variables. It is not altered by the routine. Restriction: N>0.
c 
c-----M       
c
c     INTEGER variable that must be set by the user to the number of 
c     corrections used in the BFGS update. It is not altered by the 
c     routine. Values of M less than 3 are not recommended; large 
c     values of M will result in excessive computing time. 
c     3<= M <=7 is recommended. Restriction: M>0.
c 
c-----X(1..N)
c
c     DOUBLE PRECISION array of length N. It must be set by the user 
c     to the values of the initial estimate of the solution vector. 
c     On exit with iflag=0, it contains the values of the variables 
c     at the best point found (usually a solution).
c 
c-----DFU
c
c     Name of the routine that computes the function that is minimized
c     and its gradient.
c 
c-----QPR
c
c     Variable that defines the quantity of output.
c     Depending on the value of QPR the array IPRINT(1..2) is
c     initializated in a different form.
c
c     QPR .lt. 0  --> iprint(1) = -1  iprint(2) = -1 
c     QPR .eq. 0  --> iprint(1) = -1  iprint(2) =  0
c     QPR .eq. 1  --> iprint(1) =  0  iprint(2) =  0
c     QPR .eq. 2  --> iprint(1) =  0  iprint(2) =  1
c     QPR .eq. 3  --> iprint(1) =  1  iprint(2) =  2
c     QPR .ge. 4  --> iprint(1) =  2  iprint(2) =  3
c
c     IPRINT(1) specifies the frequency of the output:
c      < 0 : no output is generated,
c      = 0 : output only at first and last iteration,
c      > 0 : output every IPRINT(1) iterations.
c 
c     IPRINT(2) specifies the type of output generated:
c      = -1 : Nothing is written. 
c      = 0 : iteration count, number of function evaluations,
c            function value, norm of the gradient, and steplength,
c      = 1 : same as IPRINT(2)=0, plus vector of variables and gradient
c            vector at the initial point,
c      = 2 : same as IPRINT(2)=1, plus vector of variables,
c      = 3 : same as IPRINT(2)=2, plus gradient vector.
c 
c-----EPS     
c
c     Positive DOUBLE PRECISION variable that must be set by the user, 
c     and determines the accuracy with which the solution is to be 
c     found. The subroutine terminates when |G| < EPS max(1,|X|), 
c     where |.| denotes the Euclidean norm.
c
c-----LOUT 
c
c     Logical output unit.
c
c-----LERR 
c
c     Logical unit for the printing of error messages.
c     (no error messages are written if lerr < 0).
c
c-----MAXIT
c
c     Maximum number of iterations.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     DESCRIPTION OF THE OUTPUT PARAMETERS.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     X(1..N) ---> The final set of (usually optimized) variables.
c     F       ---> The value of F(X).
c     G(1..N) ---> The gradient at X().
c     ITER    ---> Number of iterations used in the optimization.
c     warn    ---> Error code = .fals. if the run ends OK.
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      include 'implicit.inc'
      include 'error.inc'
      parameter (msave=7)
      integer iprint(2),iflag,icall,n,m,mp,lp,qpr
      real(kind=8), allocatable,dimension (:)  :: diag,w
      real(kind=8) x(n),g(n)
      real(kind=8) f,eps,xtol
      logical diagco,converge,warn
      real(kind=8) gtol,stpmin,stpmax
      external dfu
c
      call timer (2,ifdlbfgs,'_fdlbfgs  ',-1)  
      allocate (diag(n))
      allocate (w(n*(2*msave+1)+2*msave))
      mp = lout
      lp = lerr
      if (qpr.lt.0) then
        iprint(1)=-1
        iprint(2)=-1
      else if (qpr.eq.0) then
        iprint(1)=-1
        iprint(2)= 0
      else if (qpr.eq.1) then
        iprint(1)= 0
        iprint(2)= 0
      else if (qpr.eq.2) then
        iprint(1)= 0
        iprint(2)= 1
      else if (qpr.eq.3) then
        iprint(1)= 1
        iprint(2)= 2
      else
        iprint(1)= 2
        iprint(2)= 3
      endif
c
c-----Fill the data in 'machdep.inc' file.
c
      gtol   = 9D-1
      stpmin = 1d-20
      stpmax = 1d+20
c
c--------------------------------------------------------------------
c-----GTOL is a DOUBLE PRECISION variable with default value 0.9, 
c     which controls the accuracy of the line search routine MCSRCH. 
c     If the function and gradient evaluations are inexpensive with 
c     respect to the cost of the iteration (which is sometimes the 
c     case when solving very large problems) it may be advantageous 
c     to set GTOL to a small value. A typical small value is 0.1. 
c     Restriction: GTOL should be greater than 1.D-04.
c
c-----STPMIN and STPMAX are non-negative DOUBLE PRECISION variables 
c     which specify lower and uper bounds for the step in the line 
c     search. Their default values are 1.D-20 and 1.D+20, respectively.
c     These values need not be modified unless the exponents are too 
c     large for the machine being used, or unless the problem is 
c     extremely badly scaled (in which case the exponents should be 
c     increased).
c
c--------------------------------------------------------------------
c
c-----We do not wish to provide the diagonal matrices Hk0, and 
c     therefore set DIAGCO to FALSE.
c
      diagco = .false.
      dguess = 0.1D0
c
c-----Machine precission. XTOL is a positive DOUBLE PRECISION 
c     variable that must be set by the user to an estimate of the 
c     machine precision (e.g.10**(-16) on a SUN station 3/60). The 
c     line search routine will terminate if the relative width of the 
c     interval of uncertainty is less than XTOL.       
c
      xtol     = 2.220446049250313D-16
      icall    = 1
      iflag    = 0
      converge = .false.
      warn     = .false.
      do while (icall.le.maxit .and. (.not.converge))
        call dfu (n,x,f,g)
c
c-------Convergence test on the modulus of gradient. !!! Note that this
c       test if done even if lbfgs routine has returned a value of 
c       IFLAG equal to 1
c
        call lbfgs (n,m,x,f,g,diagco,diag,iprint,dguess,eps,xtol,w,
     &    iflag,gtol,stpmin,stpmax,mp,lp)
        if (iflag.lt.0) then
          converge = .true.
        else if (iflag.eq.1) then
          icall = icall + 1
        else if (iflag.eq.0) then
          converge = .true.
          iter = icall
          warn = .false.
          deallocate (diag)
          deallocate (w)
          return
        else
          call error ('FDLBFGS','Returned unknown iflag value',faterr)
        endif
      enddo
c
      iter = icall
      warn = .false.
      if (iflag.ne.0) warn = .true.
      deallocate (diag)
      deallocate (w)
      call timer (4,ifdlbfgs,'_fdlbfgs  ',-1)  
      return
      end
