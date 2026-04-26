C
C-----------------------------------------------------------------------
C
      SUBROUTINE MYDGECO(A,LDA,N,IPVT,RCOND,Z)
      INTEGER LDA,N,IPVT(*)
      real*8 A(LDA,*),Z(*)
      real*8 RCOND
C
C     DGECO FACTORS A real*8 MATRIX BY GAUSSIAN ELIMINATION
C     AND ESTIMATES THE CONDITION OF THE MATRIX.
C
C     IF  RCOND  IS NOT NEEDED, DGEFA IS SLIGHTLY FASTER.
C     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL.
C     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI.
C     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI.
C
C     ON ENTRY
C
C        A       real*8(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        RCOND   real*8
C                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
C                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
C                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
C                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
C                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
C                           1.0 + RCOND .EQ. 1.0
C                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
C                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
C                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
C                UNDERFLOWS.
C
C        Z       real*8(N)
C                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
C                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
C                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
C                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     LINPACK DGEFA
C     BLAS MYDAXPY,DDOT,DSCAL,DASUM
C     FORTRAN DABS,DMAX1,DSIGN
C
C     INTERNAL VARIABLES
C
      real*8 DDOT,EK,T,WK,WKM
      real*8 ANORM,S,DASUM,SM,YNORM
      INTEGER INFO,J,K,KB,KP1,L
C
C
C     COMPUTE 1-NORM OF A
C
      ANORM = 0.0D0
      DO 10 J = 1, N
         ANORM = DMAX1(ANORM,DASUM(N,A(1,J),1))
   10 CONTINUE
C
C     FACTOR
C
      CALL DGEFA(A,LDA,N,IPVT,INFO)
C
C     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
C     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
C     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
C     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
C     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
C     OVERFLOW.
C
C     SOLVE TRANS(U)*W = E
C
      EK = 1.0D0
      DO 20 J = 1, N
         Z(J) = 0.0D0
   20 CONTINUE
      DO 100 K = 1, N
         IF (Z(K) .NE. 0.0D0) EK = DSIGN(EK,-Z(K))
         IF (DABS(EK-Z(K)) .LE. DABS(A(K,K))) GO TO 30
            S = DABS(A(K,K))/DABS(EK-Z(K))
            CALL DSCAL(N,S,Z,1)
            EK = S*EK
   30    CONTINUE
         WK = EK - Z(K)
         WKM = -EK - Z(K)
         S = DABS(WK)
         SM = DABS(WKM)
         IF (A(K,K) .EQ. 0.0D0) GO TO 40
            WK = WK/A(K,K)
            WKM = WKM/A(K,K)
         GO TO 50
   40    CONTINUE
            WK = 1.0D0
            WKM = 1.0D0
   50    CONTINUE
         KP1 = K + 1
         IF (KP1 .GT. N) GO TO 90
            DO 60 J = KP1, N
               SM = SM + DABS(Z(J)+WKM*A(K,J))
               Z(J) = Z(J) + WK*A(K,J)
               S = S + DABS(Z(J))
   60       CONTINUE
            IF (S .GE. SM) GO TO 80
               T = WKM - WK
               WK = WKM
               DO 70 J = KP1, N
                  Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
         Z(K) = WK
  100 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
C     SOLVE TRANS(L)*Y = W
C
      DO 120 KB = 1, N
         K = N + 1 - KB
         IF (K .LT. N) Z(K) = Z(K) + DDOT(N-K,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 110
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
  110    CONTINUE
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
  120 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
C
      YNORM = 1.0D0
C
C     SOLVE L*V = Y
C
      DO 140 K = 1, N
         L = IPVT(K)
         T = Z(L)
         Z(L) = Z(K)
         Z(K) = T
         IF (K .LT. N) CALL MYDAXPY(N-K,T,A(K+1,K),1,Z(K+1),1)
         IF (DABS(Z(K)) .LE. 1.0D0) GO TO 130
            S = 1.0D0/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  130    CONTINUE
  140 CONTINUE
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
C     SOLVE  U*Z = V
C
      DO 160 KB = 1, N
         K = N + 1 - KB
         IF (DABS(Z(K)) .LE. DABS(A(K,K))) GO TO 150
            S = DABS(A(K,K))/DABS(Z(K))
            CALL DSCAL(N,S,Z,1)
            YNORM = S*YNORM
  150    CONTINUE
         IF (A(K,K) .NE. 0.0D0) Z(K) = Z(K)/A(K,K)
         IF (A(K,K) .EQ. 0.0D0) Z(K) = 1.0D0
         T = -Z(K)
         CALL MYDAXPY(K-1,T,A(1,K),1,Z(1),1)
  160 CONTINUE
C     MAKE ZNORM = 1.0
      S = 1.0D0/DASUM(N,Z,1)
      CALL DSCAL(N,S,Z,1)
      YNORM = S*YNORM
C
      IF (ANORM .NE. 0.0D0) RCOND = YNORM/ANORM
      IF (ANORM .EQ. 0.0D0) RCOND = 0.0D0
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MYDGEDI(A,LDA,N,IPVT,DET,WORK,JOB)
      INTEGER LDA,N,IPVT(*),JOB
      real*8 A(LDA,*),DET(2),WORK(*)
C
C     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX
C     USING THE FACTORS COMPUTED BY DGECO OR DGEFA.
C
C     ON ENTRY
C
C        A       real*8(LDA, N)
C                THE OUTPUT FROM DGECO OR DGEFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR DGEFA.
C
C        WORK    real*8(N)
C                WORK VECTOR.  CONTENTS DESTROYED.
C
C        JOB     INTEGER
C                = 11   BOTH DETERMINANT AND INVERSE.
C                = 01   INVERSE ONLY.
C                = 10   DETERMINANT ONLY.
C
C     ON RETURN
C
C        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE UNCHANGED.
C
C        DET     real*8(2)
C                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED.
C                OTHERWISE NOT REFERENCED.
C                DETERMINANT = DET(1) * 10.0**DET(2)
C                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0
C                OR  DET(1) .EQ. 0.0 .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS
C        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED.
C        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY
C        AND IF DGECO HAS SET RCOND .GT. 0.0 OR DGEFA HAS SET
C        INFO .EQ. 0 .
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS DAXPY,DSCAL,DSWAP
C     FORTRAN DABS,MOD
C
C     INTERNAL VARIABLES
C
      real*8 T
      real*8 TEN
      INTEGER I,J,K,KB,KP1,L,NM1
C
C
C     COMPUTE DETERMINANT
C
      IF (JOB/10 .EQ. 0) GO TO 70
         DET(1) = 1.0D0
         DET(2) = 0.0D0
         TEN = 10.0D0
         DO 50 I = 1, N
            IF (IPVT(I) .NE. I) DET(1) = -DET(1)
            DET(1) = A(I,I)*DET(1)
C        ...EXIT
            IF (DET(1) .EQ. 0.0D0) GO TO 60
   10       IF (DABS(DET(1)) .GE. 1.0D0) GO TO 20
               DET(1) = TEN*DET(1)
               DET(2) = DET(2) - 1.0D0
            GO TO 10
   20       CONTINUE
   30       IF (DABS(DET(1)) .LT. TEN) GO TO 40
               DET(1) = DET(1)/TEN
               DET(2) = DET(2) + 1.0D0
            GO TO 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     COMPUTE INVERSE(U)
C
      IF (MOD(JOB,10) .EQ. 0) GO TO 150
         DO 100 K = 1, N
            A(K,K) = 1.0D0/A(K,K)
            T = -A(K,K)
            CALL DSCAL(K-1,T,A(1,K),1)
            KP1 = K + 1
            IF (N .LT. KP1) GO TO 90
            DO 80 J = KP1, N
               T = A(K,J)
               A(K,J) = 0.0D0
               CALL MYDAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
C
C        FORM INVERSE(U)*INVERSE(L)
C
         NM1 = N - 1
         IF (NM1 .LT. 1) GO TO 140
         DO 130 KB = 1, NM1
            K = N - KB
            KP1 = K + 1
            DO 110 I = KP1, N
               WORK(I) = A(I,K)
               A(I,K) = 0.0D0
  110       CONTINUE
            DO 120 J = KP1, N
               T = WORK(J)
               CALL MYDAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
            L = IPVT(K)
            IF (L .NE. K) CALL DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      real*8 FUNCTION DASUM(N,DX,INCX)
C
C     TAKES THE SUM OF THE ABSOLUTE VALUES.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      real*8 DX(*),DTEMP
      INTEGER I,INCX,M,MP1,N,NINCX
C
      DASUM = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX = N*INCX
      DO 10 I = 1,NINCX,INCX
        DTEMP = DTEMP + DABS(DX(I))
   10 CONTINUE
      DASUM = DTEMP
      RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,6)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DABS(DX(I))
   30 CONTINUE
      IF( N .LT. 6 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,6
        DTEMP = DTEMP + DABS(DX(I)) + DABS(DX(I + 1)) + DABS(DX(I + 2))
     *  + DABS(DX(I + 3)) + DABS(DX(I + 4)) + DABS(DX(I + 5))
   50 CONTINUE
   60 DASUM = DTEMP
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DGEFA(A,LDA,N,IPVT,INFO)
      INTEGER LDA,N,IPVT(*),INFO
      real*8 A(LDA,*)
C
C     DGEFA FACTORS A real*8 MATRIX BY GAUSSIAN ELIMINATION.
C
C     DGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR DGEFA) .
C
C     ON ENTRY
C
C        A       real*8(LDA, N)
C                THE MATRIX TO BE FACTORED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
C                WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C                = 0  NORMAL VALUE.
C                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
C                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO
C                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS MYDAXPY,DSCAL,IDAMAX
C
C     INTERNAL VARIABLES
C
      real*8 T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL MYDAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
C
C-----------------------------------------------------------------------
C
*     SUBROUTINE  DSCAL(N,DA,DX,INCX)
C
C     SCALES A VECTOR BY A CONSTANT.
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
*     real*8 DA,DX(*)
*     INTEGER I,INCX,M,MP1,N,NINCX
C
*     IF(N.LE.0)RETURN
*     IF(INCX.EQ.1)GO TO 20
C
C        CODE FOR INCREMENT NOT EQUAL TO 1
C
*     NINCX = N*INCX
*     DO 10 I = 1,NINCX,INCX
*       DX(I) = DA*DX(I)
*  10 CONTINUE
*     RETURN
C
C        CODE FOR INCREMENT EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
*  20 M = MOD(N,5)
*     IF( M .EQ. 0 ) GO TO 40
*     DO 30 I = 1,M
*       DX(I) = DA*DX(I)
*  30 CONTINUE
*     IF( N .LT. 5 ) RETURN
*  40 MP1 = M + 1
*     DO 50 I = MP1,N,5
*       DX(I) = DA*DX(I)
*       DX(I + 1) = DA*DX(I + 1)
*       DX(I + 2) = DA*DX(I + 2)
*       DX(I + 3) = DA*DX(I + 3)
*       DX(I + 4) = DA*DX(I + 4)
*  50 CONTINUE
*     RETURN

C      real*8 FUNCTION DDOT(N,DX,INCX,DY,INCY)
CC
CC     FORMS THE DOT PRODUCT OF TWO VECTORS.
CC     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
CC     JACK DONGARRA, LINPACK, 3/11/78.
CC
C      real*8 DX(*),DY(*),DTEMP
C      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
CC
C      DDOT = 0.0D0
C      DTEMP = 0.0D0
C      IF(N.LE.0)RETURN
C      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
CC
CC        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
CC          NOT EQUAL TO 1
CC
C      IX = 1
C      IY = 1
C      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
C      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
C      DO 10 I = 1,N
C        DTEMP = DTEMP + DX(IX)*DY(IY)
C        IX = IX + INCX
C        IY = IY + INCY
C   10 CONTINUE
C      DDOT = DTEMP
C      RETURN
CC
CC        CODE FOR BOTH INCREMENTS EQUAL TO 1
CC
CC
CC        CLEAN-UP LOOP
CC
C   20 M = MOD(N,5)
C      IF( M .EQ. 0 ) GO TO 40
C      DO 30 I = 1,M
C        DTEMP = DTEMP + DX(I)*DY(I)
C   30 CONTINUE
C      IF( N .LT. 5 ) GO TO 60
C   40 MP1 = M + 1
C      DO 50 I = MP1,N,5
C        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
C     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
C   50 CONTINUE
C   60 DDOT = DTEMP
C      RETURN
C      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MYDAXPY(N,DA,DX,INCX,DY,INCY)
      real*8 DX(*),DY(*),DA
      INTEGER I,INCX,INCY,M,MP1,N
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE  DSWAP (N,DX,INCX,DY,INCY)
C
C     INTERCHANGES TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      real*8 DX(*),DY(*),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
C         TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DX(IX)
        DX(IX) = DY(IY)
        DY(IY) = DTEMP
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C       CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C       CLEAN-UP LOOP
C
   20 M = MOD(N,3)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
   30 CONTINUE
      IF( N .LT. 3 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,3
        DTEMP = DX(I)
        DX(I) = DY(I)
        DY(I) = DTEMP
        DTEMP = DX(I + 1)
        DX(I + 1) = DY(I + 1)
        DY(I + 1) = DTEMP
        DTEMP = DX(I + 2)
        DX(I + 2) = DY(I + 2)
        DY(I + 2) = DTEMP
   50 CONTINUE
      RETURN
      END
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
      real*8 DX(*),DMAX
      INTEGER I,INCX,IX,N
      IDAMAX = 0
      IF( N .LT. 1 ) RETURN
      IDAMAX = 1
      IF(N.EQ.1)RETURN
      IF(INCX.EQ.1)GO TO 20
      IX = 1
      DMAX = DABS(DX(1))
      IX = IX + INCX
      DO 10 I = 2,N
         IF(DABS(DX(IX)).LE.DMAX) GO TO 5
         IDAMAX = I
         DMAX = DABS(DX(IX))
    5    IX = IX + INCX
   10 CONTINUE
      RETURN
   20 DMAX = DABS(DX(1))
      DO 30 I = 2,N
         IF(DABS(DX(I)).LE.DMAX) GO TO 30
         IDAMAX = I
         DMAX = DABS(DX(I))
   30 CONTINUE
      RETURN
      END
c------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C
        SUBROUTINE      TRACE (H, E, W, N, IERR)
C
C TRACE CALLS TREDIG AND TLQGRM TO DIAGONALIZE A SYMMETRIC REAL MATRIX.
C THE MATRIX IS PASSED DOWN IN H AND IS REPLACED BY THE EIGENVECTORS.
C THE EIGENVALUES IN E ARE STORED SMALLEST FIRST.
C THE WORK STORE W SHOULD BE AT LEAST OF DIMENSION N.
C SKK ==================================================================

        IMPLICIT        real*8 (A-H, O-Z)
        DIMENSION       H(N,N), E(*), W(*)

        CALL TREDIG     (N, E, W, H)
        CALL TQLGRM     (N, E, W, H, IERR)

        RETURN
        END
c------------------------------------------------------------------
        SUBROUTINE      TREDIG (N, D, E, Z)
        IMPLICIT        real*8 (A-H, O-Z)
        DIMENSION       D(*), E(*), Z(N,N)
        PARAMETER       (ZERO = 0.0D0, ONE = 1.0D0)

        IF (N .EQ. 1) GOTO 320

        DO 30 II = 2,N
        I       = N + 2 - II
        L       = I - 1
        H       = ZERO
        SCALE   = ZERO

        IF (L .LT. 2) GOTO 130

                DO 31 K = 1,L
                SCALE   = SCALE + DABS(Z(I,K))
31              CONTINUE

        IF (SCALE .NE. ZERO) GOTO 140
130     E(I)    = Z(I,L)
        GOTO 290

140     RSCALE = ONE/SCALE
                DO 32 K = 1,L
                Z(I,K)  = Z(I,K)*RSCALE
                H       = H + Z(I,K)*Z(I,K)
32              CONTINUE
        F       = Z(I,L)
        G       = -DSIGN(DSQRT(H),F)
        E(I)    = SCALE*G
        H       = H - F*G
        Z(I,L)  = F - G
        F       = ZERO
        RH      = ONE/H
        RHSCALE = RH*RSCALE        

                DO 33 J = 1,L
                Z(J,I)  = Z(I,J)*RHSCALE
                G       = ZERO

                        DO 34 K = 1,J
                        G       = G + Z(J,K)*Z(I,K)
34                      CONTINUE

                JP1 = J + 1
                IF (L .LT. JP1) GOTO 220

                        DO 35 K = JP1,L
                        G       = G + Z(K,J)*Z(I,K)
35                      CONTINUE

220             E(J)    = G*RH
                F       = F + E(J)*Z(I,J)
33              CONTINUE

        HH = F/(H + H)

                DO 36 J = 1,L
                F       = Z(I,J)
                G       = E(J) - HH*F
                E(J)    = G
                        DO 37 K = 1,J
                        Z(J,K)  = Z(J,K) - F*E(K) - G*Z(I,K)
37                      CONTINUE
36              CONTINUE

                DO 38 K        = 1,L
                Z(I,K)  =  SCALE*Z(I,K)
38              CONTINUE

290     D(I)    = H
30      CONTINUE

320     D(1)    = ZERO
        E(1)    = ZERO

        DO 500 I = 1,N
        L       = I - 1
        IF (D(I) .EQ. ZERO) GOTO 380

                DO 40 J = 1,L
                G       = ZERO

                        DO 41 K = 1,L
                        G       = G + Z(I,K)*Z(K,J)
41                      CONTINUE

                        DO 42 K = 1,L
                        Z(K,J)  = Z(K,J) - G*Z(K,I)
42                      CONTINUE

40              CONTINUE

380     D(I)    = Z(I,I)
        Z(I,I)  = ONE
        IF(L .LT. 1) GOTO 500

                DO 43 J = 1,L
                Z(J,I)  = ZERO
                Z(I,J)  = ZERO
43              CONTINUE

500     CONTINUE
        RETURN
        END
C
C-----------------------------------------------------------------------
C
        SUBROUTINE      TQLGRM        (N, D, E, Z, IERR)
        IMPLICIT        real*8 (A-H, O-Z)
        DIMENSION       D(*), E(*), Z(N,N)
        PARAMETER (AMACH = 16.0E-13)
        PARAMETER (ZERO = 0.0D0, ONE = 1.0D0)
C
        IERR    = 0
        IF (N .EQ. 1) RETURN
C
        DO 30 I = 2,N
        E(I-1)  = E(I)
30      CONTINUE

        F       = ZERO
        B       = ZERO
        E(N)    = ZERO

        DO 31 L = 1,N
        J       = 0
        H       = AMACH*(DABS(D(L)) + DABS(E(L)))
        IF (B .LT. H) B = H

                DO 32 M = L,N
                IF (DABS(E(M)) .LE. B) GOTO 120
32              CONTINUE

120     IF (M .EQ. L) GOTO 220

130     IF (J .EQ. 30) THEN
                IERR    = L
                RETURN
        END IF

        J       = J + 1
        L1      = L + 1
        G       = D(L)
        P       = (D(L1) - G)/(2*E(L))
        IF (DABS(P*AMACH) .GT. ONE) THEN
                R       = P
        ELSE
                R       = DSQRT(P*P + 1)
        END IF
        D(L)    = E(L)/(P + DSIGN(R,P))
        H       = G - D(L)

                DO 33 I = L1,N
                D(I)    = D(I) - H
33              CONTINUE

        F       = F + H
        P       = D(M)
        C       = ONE
        S       = ZERO
        MML     = M - L

                DO 34 II = 1,MML
                I       = M - II
                G       = C*E(I)
                H       = C*P
                IF (DABS(P) .GE. DABS(E(I))) THEN
                        C       = E(I)/P
                        R       = DSQRT(C*C + 1)
                        E(I+1)  = S*P*R
                        S       = C/R
                        C       = ONE/R
                ELSE
                        C       = P/E(I)
                        R       = DSQRT(C*C + 1)
                        E(I+1)  = S*E(I)*R
                        S       = 1.D0/R
                        C       = C*S
                END IF
                P       = C*D(I) - S*G
                D(I+1)  = H + S*(C*G + S*D(I))

                        DO 35 K = 1,N
                        H       = Z(K,I+1)
                        Z(K,I+1)= S*Z(K,I) + C*H
                        Z(K,I)  = C*Z(K,I) - S*H
35                      CONTINUE

34              CONTINUE

        E(L)    = S*P
        D(L)    = C*P
        IF (DABS(E(L)) .GT. B) GOTO 130

220     D(L)    = D(L) + F
31      CONTINUE

        DO 300 II = 2,N
        I       = II - 1
        K       = I
        P       = D(I)

                DO 260 J = II,N
                IF (D(J) .GE. P) GOTO 260
                K       = J
                P       = D(J)
260             CONTINUE

        IF (K .EQ. I) GOTO 300
        D(K)    = D(I)
        D(I)    = P

                DO 37 J = 1,N
                P       = Z(J,I)
                Z(J,I)  = Z(J,K)
                Z(J,K)  = P
37              CONTINUE

300     CONTINUE
        RETURN
        END


