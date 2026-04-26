C-----------------------------------------------------------------------
C
C     FINDSYMMETRY: A program designed to locate symmetry elements and 
C                   to symmetrize molecules.
C     Authors: Beruski, O. and Vidal, L.N.
C     Contact: lnvidal@utfpr.edu.br
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NEWSYM (NCENT,LUW,NEQ,XYZMOL,CMOL,INEQM,MULTM,ATNAM)
C
      IMPLICIT REAL(KIND=8) (A-H,O-Z)
      INTEGER(KIND=4)     NCENT
      CHARACTER           LINE*80, PG*4, ATOMNAME*2, ANSWER*1
      CHARACTER(LEN=24)   FFDATE
      CHARACTER(LEN=40)   HSTNAM
C
C     BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccc
C
      REAL(KIND=8) XYZMOL(NCENT,3)
      REAL(KIND=8) CMOL(NCENT)
      CHARACTER(LEN=8) ATNAM(NCENT)

      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR

      INTEGER   INEQ(MAXAT),MULT(MAXAT)
      INTEGER   INEQM(NCENT),MULTM(NCENT)
C
C-----------------------------------------------------------------------
C
      DIMENSION XYZorig(MAXAT,3)

      IF (NCENT.GT.MAXAT) THEN
        WRITE (0,*)'# NCENT = ',NCENT,' > MAXAT = ',MAXAT
        STOP ' # newsym.f: Increase the value of MAXAT'
      ENDIF
      WRITE (LUW,111)
 111  FORMAT (' #',/,
     !' # ---------------------------------------------',/,' #',/,
     !' # B E G I N   S Y M M E T R Y   A N A L Y S I S',/,' #',/,
     !' # ---------------------------------------------',/,' #')

      AU = .TRUE.
      NATOM = NCENT
      XYZ(1:NATOM,1:3) = XYZMOL(1:NATOM,1:3)
      CHARGE(1:NATOM)  = CMOL(1:NATOM)
      DO I=1,NATOM
        ATMSBL(I) = ATNAM(I)(3:4)
        ISOTOPE   = 1
        RMASS(I)  = ABISTP(NINT(CHARGE(I)),ISOTOPE,"MASS")
      ENDDO
      XYZorig = XYZ
      CALL ROTOANA (LUW)
      MAXCYCL =  2
      answer  = "n"

 01   CONTINUE

      XYZ = XYZorig
      CALL XYZzero
      IF ( answer.EQ."y" ) CALL SYMMETRIZER (MAXCYCL)
      CALL FINDPG (LUW,NCENT,NEQ,INEQ,MULT,PG)
      INEQM(1:NCENT) = INEQ(1:NCENT)
      MULTM(1:NCENT) = MULT(1:NCENT)
C
C    "SYMMETRIZER" MODULE.
C
      IF ( PG.EQ."C1") THEN
         MAXCYCL=MAXCYCL-1
         answer = "y"
         WRITE (*,*) "PG is C1 --> Trying to symmetrize the molecule"
      END IF 
      IF ( answer.EQ."y" .AND. MAXCYCL.GT.0 ) GO TO 01
C
C     Date and Time
C
      CALL FDATE  (FFDATE)
      CALL HOSTNM (HSTNAM)
      WRITE (LUW,'(/1X,2A/1X,2A/)')
     !'Host name     : ',HSTNAM,
     !'Date and time : ',FFDATE

*     WRITE (LUW,'(A/A/A/)')
*    !' ***********************',
*    !' * END OF CALCULATIONS *',
*    !' ***********************'

*     CLOSE (80)
      WRITE (LUW,112)
 112  FORMAT (' #',/,
     !' # -------------------------------------------',/,' #',/,
     !' # D O N E   S Y M M E T R Y   A N A L Y S I S',/,' #',/,
     !' # -------------------------------------------',/,' #')

 05   CONTINUE

*     WRITE (*,"(//A/A//)") 
*    *" The program found an error while reading atomic coordinates.",
*    *" Check if the number of atoms is correct. "

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Module for automatic detection of molecular point group
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FINDPG(LUW,NCENT,NEQ,INEQ,MULT,PG)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION DAiAj(MAXAT,MAXAT), Di(MAXAT), Dj(MAXAT),
     *          vec(3), SMTZ(3,3), SIGMAVECS(MAXOPSYM,3), 
     *          SYMVEC(MAXOPSYM,3), VA(3),VB(3),VC(3),
     *          C2PRP(MAXOPSYM,5),CNAXIS(MAXOPSYM,5),
     *          SNAXIS(MAXOPSYM,5)

      INTEGER   INEQ(MAXAT),MULT(MAXAT)

      CHARACTER PG*4
      LOGICAL   SEA(MAXAT), CONFIRM, INVERSION, SIGMAV, SPHERICAL

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0

      SPHERICAL=.FALSE.

      IF ( NATOM.EQ.1 ) THEN
         WRITE (LUW,'(/A/)') " * Point Group found : Kh "
         RETURN
      END IF

*     WRITE (LUW,'(/A/A/A/)')
*    !' *****************************************',
*    !' **** OUTPUT FROM PROGRAM POINT GROUP ****',
*    !' *****************************************'

c Threshold to consider two distances equal
      Thres = 5.D-3
      IF (AU) Thres = Thres / BOHR
      WRITE (LUW,'(A,1P,D8.2/)')
     !" * Threshold to consider two distances equal (angstrom) : ",Thres

c Threshold to consider two inertia moments equal
      ThresIABC = 1.D-1
      WRITE (LUW,'(A,1P,D8.2/)')
     !" * Threshold to consider two inertia moments equal : ",ThresIABC

c Center of Mass translation
      CMX = 0.D0
      CMY = 0.D0
      CMZ = 0.D0
      TMASS = 0.D0
      DO III = 1,NATOM
         AMASS = RMASS(III)
         CMX = CMX + XYZ(III,1) * AMASS
         CMY = CMY + XYZ(III,2) * AMASS
         CMZ = CMZ + XYZ(III,3) * AMASS
         TMASS = TMASS + AMASS
      END DO
      CMX = CMX / TMASS
      CMY = CMY / TMASS
      CMZ = CMZ / TMASS

c printing atomic coordinates
      WRITE (LUW,'(/A,A/A/)')
     !' * The following coordinates were used to locate',
     !' the point group : ',
     !'   (center of mass in the origin)'

      CONV = 1.D0
      IF ( AU ) CONV = BOHR
      WRITE (LUW,'(A,A,13X,A,13X,A)')
     !'   #  Atom  Z   Isotopic Mass        ','x','y','z'
      DO III = 1,NATOM
         AMASS = RMASS(III)
         XYZ(III,1) = ( XYZ(III,1)-CMX ) * CONV
         XYZ(III,2) = ( XYZ(III,2)-CMY ) * CONV
         XYZ(III,3) = ( XYZ(III,3)-CMZ ) * CONV
         WRITE (LUW,'(1X,I3,2X,A4,F4.0,3X,F12.6,2X,3F14.6)')
     !III,ATMSBL(III),CHARGE(III),AMASS,(XYZ(III,J),J=1,3)
      END DO

c Interatomic distances (in Angstrom)
      DO I = 1,NATOM
         DO J = I,NATOM
            XX = ( XYZ(I,1) - XYZ(J,1) ) * ( XYZ(I,1) - XYZ(J,1) )
            YY = ( XYZ(I,2) - XYZ(J,2) ) * ( XYZ(I,2) - XYZ(J,2) )
            ZZ = ( XYZ(I,3) - XYZ(J,3) ) * ( XYZ(I,3) - XYZ(J,3) )
            DAiAj(I,J) = DSQRT( XX + YY + ZZ )
            DAiAj(J,I) = DAiAj(I,J)
            IF ( I.NE.J .AND. DAiAj(I,J).LT.Thres ) THEN
               WRITE (LUW,"(//A,I3,A,I3,A//)")
     *" Error: Atoms ",I," and ",J," are too close! "
               STOP
            END IF
         END DO
      END DO
      WRITE (LUW,'(/A/A)')
     !' Interatomic distances (in A) ',
     !' ---------------------------- '
c printing up to 14 columns
      MAXCOL = 8
      IF ( NATOM.LT.MAXCOL ) MAXCOL = NATOM
      NCOLi  = 1
      NCOL   = MAXCOL
      DO WHILE ( NCOL.LE.NATOM )
         WRITE (LUW,'(10X,14(I3,1X,A2,3X))') 
     !   (i,ATMSBL(i),i=NCOLi,NCOL)
         DO i = 1,NATOM
            WRITE (LUW,'(1X,I3,1X,A2,1X,14(F8.4,1X))') 
     !i,ATMSBL(i),(DAiAj(i,j),j=NCOLi,NCOL)
         END DO
         NCOL  = NCOL  + MAXCOL
         NCOLi = NCOLi + MAXCOL
         WRITE (LUW,*)
      END DO
      NCOL  = NCOL  - MAXCOL
      NCOLi = NCOLi - MAXCOL
      IF ( NCOL.LT.NATOM ) THEN
         NCOLi = NCOL + 1
         NCOL = NATOM
         WRITE (LUW,'(10X,14(I3,1X,A2,3X))') 
     !   (i,ATMSBL(i),i=NCOLi,NCOL)
         DO i = 1,NATOM
            WRITE (LUW,'(1X,I3,1X,A2,1X,14(F8.4,1X))') 
     !i,ATMSBL(i),(DAiAj(i,j),j=NCOLi,NCOL)
         END DO
      END IF

c Separating atoms into sets of Symmetrically Equivalent Ones (SEA)
c      ISEA( set , atoms into set, e.g., atoms 1,4,12)
c      NSEA( set ) = number of atoms in the set
c      SEA ( natom ) : if .TRUE.  , atom "i" is symmetrically-equivalent
c                      if .FALSE. , atom "i" is symmetry-unique 
      DO i=1,MAXAT
         DO j=1,MAXAT
            ISEA(i,j) = 0
         END DO
         SEA(i)    = .FALSE.
      END DO
      nseaset = 0
      DO i = 1,NATOM
         IF ( .NOT.SEA(i) ) THEN
            nseaset = nseaset + 1
            NSEA(nseaset) = 1
            ISEA(nseaset,NSEA(nseaset)) = i
            SEA(i) = .TRUE.
         END IF
         DO ii=1,NATOM
            Di(ii) = DAiAj(i,ii)
         END DO
         j=i+1
         DO WHILE ( j.LE.NATOM ) 
            NEQD = 0
            DO jj=1,NATOM
               Dj(jj) = DAiAj(j,jj)
            END DO
            RNORMDIFF = DABS( DiNORM(Di,NATOM,MAXAT) -
     &                        DiNORM(Dj,NATOM,MAXAT) )
            IF ( RNORMDIFF.LE.Thres .AND. RMASS(i).EQ.RMASS(j) .AND.
     &           .NOT.SEA(j) ) THEN
               DO ii=1,NATOM
               DO jj=1,NATOM
                  IF ( DABS(Di(ii)-Dj(jj)).LE.Thres  .AND.
     &                 RMASS(ii).EQ.RMASS(jj)) THEN
                      NEQD = NEQD+1 
                      Dj(jj) = -1.D6
                  END IF
               END DO
               END DO
            END IF
            IF ( NEQD.EQ.NATOM ) THEN
               NSEA(nseaset) = NSEA(nseaset) + 1
               ISEA(nseaset,NSEA(nseaset)) = j
               SEA(j) = .TRUE.
            END IF
            j=j+1
         END DO
      END DO

      MAXCOL=15
      WRITE (LUW,'(/A/A/A/A)') 
     !" Symmetrically Equivalent Atoms (SEA): ",
     !" -------------------------------------   ",
     !" Set | # Atoms | SEA ",
     !" -------------------------------------   "
      DO i=1,nseaset
         Ji=1
         Jf=NSEA(i)
         IF ( Jf.LT.MAXCOL ) 
     !WRITE (LUW,'(1X,I3,1X,"|",3X,I3,3X,"|",1X,15I4)')
     !i, NSEA(i),(ISEA(i,j),j=Ji,Jf)
         IF ( Jf.GE.MAXCOL ) THEN
            Jf=MAXCOL
            DO WHILE ( Jf.LE.NSEA(i) )
               IF (Ji.EQ.1) THEN
      WRITE (LUW,'(1X,I3,1X,"|",3X,I3,3X,"|",1X,15I4)')
     !i, NSEA(i),(ISEA(i,j),j=Ji,Jf)
               ELSE
      WRITE (LUW,'(1X,3X,1X,1X,3X,3X,3X,"|",1X,15I4)') 
     !(ISEA(i,j),j=Ji,Jf)
               END IF
               Ji=Ji+MAXCOL
               Jf=Jf+MAXCOL
            END DO
            Jf=Jf-MAXCOL
         END IF
         IF ( Jf.LT.NSEA(i)) THEN
            Ji=Jf+1
            Jf=NSEA(i)
      WRITE (LUW,'(1X,3X,1X,1X,3X,3X,3X,"|",1X,15I4)') 
     !(ISEA(i,j),j=Ji,Jf)
         END IF
*        WRITE (LUW,*)
      END DO

      neq = nseaset
      do i = 1,neq
        ineq(i) = isea(i,1)
        mult(i) = nsea(i)
      enddo
*EFM
*     write (luw,82) neq, (ineq(i),i=1,neq)
*     write (luw,84) (mult(i),i=1,neq)
*82   format (1x,'# Number of Inequivalent atoms: ',i4,/
*    &        1x,'# Ineqs                    :', 1000(/,' #',40i4))
*84   format (1x,'# Multiplicity of each Ineq:',1000(/,' #',40i4))

ccccccccccccccccccccccccccccccccccccccccccccccc
c Checking if molecule has a center of symmetry
ccccccccccccccccccccccccccccccccccccccccccccccc
      INVERSION=.FALSE.
      DO i=1,3
      DO j=1,3
         SMTZ(i,j) = 0.D0
         IF (i.eq.j) SMTZ(i,j) = -1.D0
      END DO
      END DO
      CALL CHECKSYMOPER (SMTZ,Thres,CONFIRM)
      IF ( CONFIRM ) THEN
         WRITE (LUW,'(/A/A)')
     !" * This molecule possesses a center of symmetry ",
     !"   -------------------------------------------- "
         INVERSION=.TRUE.
      END IF

      CALL  MOLxIAxIBxIC(PIA,PIB,PIC,VA,VB,VC)
ccccccccccccccccccccccccccccccc
c Cheking if molecule is linear
ccccccccccccccccccccccccccccccc
      IF ( DABS(PIC-PIB).LT.ThresIABC .AND. DABS(PIA).LT.ThresIABC )THEN
         WRITE (LUW,"(/A/A/A/A,F8.5,A,F8.5,A,F8.5,A)")
     !" * Molecule recongized as linear ",
     !" --------------------------------",
     !" This system can have any Cn collinear to the axis ",
     !" (",VA(1),",",VA(2),",",VA(3),")"
         PG = "Coov"
         IF ( INVERSION ) THEN
            WRITE (LUW,"(A/)")
     !" as well as C2 and reflexion planes perpendicular to this axis "
            PG = "Dooh"
         END IF
         CALL WRTPG (LUW,PG)
         RETURN
      END IF

ccccccccccccccccccccccccccccccccccccc
c Checking if molecule is Cs (planar)
ccccccccccccccccccccccccccccccccccccc
      IF ( NSEASET.EQ.NATOM .AND. DABS(PIC-PIA-PIB).LT.ThresIABC 
     !     .AND. NATOM.GT.2 ) THEN
         WRITE (LUW,"(/A/A/A/A,F8.5,A,F8.5,A,F8.5,A)")
     !" * The molecule belongs to Cs group (planar) ",
     !" -------------------------------------------",
     !" The normal to the plane vector is ",
     !" (",VC(1),",",VC(2),",",VC(3),")"
         PG = "Cs  "
         CALL WRTPG (LUW,PG)
         RETURN
      END IF

ccccccccccccccccccccccc
c Proper Rotation Axes
ccccccccccccccccccccccc
      NSYMVEC=0
      NCNAXIS=0
      NC2=0
      NC3=0
      NC4=0
      NC5=0
      DO i=1,nseaset
         MAXC2=15
         MAXC3=20
         MAXC4=6
         MAXC5=24
ccccccccccccc
c SINGLE ATOM
ccccccccccccc
c Not used to find Cn's

cccccccccccc
c TWO ATOMS
cccccccccccc
         IF (NSEA(i).EQ.2 .AND. .NOT. LINEAR ) THEN
            NSYMVEC=NSYMVEC+1
            DO j=1,3
              SYMVEC(NSYMVEC,j) = XYZ(ISEA(i,1),j)-XYZ(ISEA(i,2),j)
            END DO
            RNSYMVEC = DSQRT( SYMVEC(NSYMVEC,1)*SYMVEC(NSYMVEC,1) +
     !                        SYMVEC(NSYMVEC,2)*SYMVEC(NSYMVEC,2) +
     !                        SYMVEC(NSYMVEC,3)*SYMVEC(NSYMVEC,3) )
            DO j=1,3
               SYMVEC(NSYMVEC,j) = SYMVEC(NSYMVEC,j) / RNSYMVEC
            END DO
c            WRITE (LUW,"(/A,I3,A/A/A/A,F8.5,A,F8.5,A,F8.5,A/)")
c     !" * The analysis of SEA set ",i," composed of 2 atoms ",
c     !" provided the following useful information: ",
c     !" A C2 in this molecule (if any) must be perpendicular to",
c     !" the vector (x,y,z) = (",
c     !SYMVEC(NSYMVEC,1),",",SYMVEC(NSYMVEC,2),",",SYMVEC(NSYMVEC,3),")"
            CALL C2AXIS (NSYMVEC,SYMVEC,I,NCNAXIS,CNAXIS,Thres)

         ELSE IF ( NSEA(i).GT.2 ) THEN 
ccccccccccccccccccccc
c THREE OR MORE ATOMS
ccccccccccccccccccccc

c Inertia tensor of the SEA set
            CALL SEAxIAxIBxIC (I,PIA,PIB,PIC,VA,VB,VC,CMX,CMY,CMZ)

ccccccccccccccccccccccccccccccccc
c Polygon case
C SEA forming a planar arrangment
ccccccccccccccccccccccccccccccccc
            IF ( DABS(PIA+PIB-PIC).LT.ThresIABC ) THEN
               N=NSEA(I)
               CALL CNCONFIRM (N,VC,Thres,NCNAXIS,CNAXIS)

ccccccccccccccccccccccccccccccccccccc
c Polyhedron case
c SEA forming a prolate symmetric top
ccccccccccccccccccccccccccccccccccccc
            ELSE IF ( DABS(PIB-PIC).LT.ThresIABC .AND. 
     !                    (PIB-PIA).GT.ThresIABC ) THEN
               N=NSEA(I)
               IF ( MOD(N,2).NE.0 ) 
     !STOP "Error in pgroup: expecting an even number but found odd one"
               N=NSEA(I) / 2
               CALL CNCONFIRM (N,VA,Thres,NCNAXIS,CNAXIS)

ccccccccccccccccccccccccccccccccccccc
c Polyhedron case
c SEA forming an oblate symmetric top
ccccccccccccccccccccccccccccccccccccc
            ELSE IF ( DABS(PIA-PIB).LT.ThresIABC .AND. 
     !                    (PIC-PIB).GT.ThresIABC ) THEN
               N=NSEA(I)
               IF ( MOD(N,2).NE.0 ) 
     !STOP "Error in pgroup: expecting an even number but found odd one"
               N=NSEA(I) / 2
               CALL CNCONFIRM (N,VC,Thres,NCNAXIS,CNAXIS)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Polyhedron case
c SEA forming an asymmetric top arrangment (orthorhombic case)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            ELSE IF ( DABS(PIA-PIB).GT.ThresIABC .AND.
     !                    (PIC-PIB).GT.ThresIABC ) THEN
              CALL CNCONFIRM (2,VA,Thres,NCNAXIS,CNAXIS)

ccccccccccccccccccccc
c Spherical Rotor
ccccccccccccccccccccc
            ELSE IF ( DABS(PIA-PIB).LT.ThresIABC .AND.
     !                    (PIC-PIB).LT.ThresIABC ) THEN
              SPHERICAL=.TRUE.
c finding C2
              CALL C2SPHER (I,NCNAXIS,CNAXIS,Thres,MAXC2,NC2)
c finding C3
              CALL  C3AXIS (I,NCNAXIS,CNAXIS,Thres,MAXC3,NC3)
c finding C4
              IF ( NC2.EQ.9 )
     *        CALL  C4AXIS (I,NCNAXIS,CNAXIS,Thres,MAXC4,NC4)
c finding C5
              IF ( NC2.EQ.15 )
     *        CALL  C5AXIS (I,NCNAXIS,CNAXIS,Thres,MAXC5,NC5)

ccccccccccccccccccccccccccccccccccc
c Is there any other possibility ?
ccccccccccccccccccccccccccccccccccc
            ELSE 
               STOP " Which case is it? "
            END IF

         END IF
      END DO

c Printing Cn's found
      IF ( NCNAXIS.GT.0 ) THEN
         WRITE (LUW,'(/A,I3/A/A)')
     !" * The number of C_n^k axes found is :",NCNAXIS,
     !"   --------------------------------------",
     !" Cn-Axis   n   k   theta         Direction (x y z)"
         DO i=1,NCNAXIS
            N = NINT(CNAXIS(i,4))
            K = NINT(CNAXIS(i,5))
            THETA = 360.D0*CNAXIS(i,5)/CNAXIS(i,4)
            WRITE (LUW,'(1X,I4,3X,2I4,4X,I4,1X,3F10.5)')
     !i,N,K,NINT(theta),(CNAXIS(i,j),j=1,3)
         END DO
      END IF

ccccccccccccccccccccccccccccccccccccc
c Searching for Perpendicular-C2 axes
ccccccccccccccccccccccccccccccccccccc
      NC2PRP=0
      IF ( NCNAXIS.GT.0 .AND. .NOT.SPHERICAL ) 
     *   CALL C2PERP (NC2PRP,C2PRP,NCNAXIS,CNAXIS,Thres,NSYMVEC,SYMVEC,
     *               IPRPCN)
c printing C2perp found
      IF ( NC2PRP.GT.0 ) THEN
         WRITE (LUW,'(/A,I3/A/A)')
     !" * The number of Perpendicular C2 axes found is :",NC2PRP,
     !"   -------------------------------------------------",
     !" C2-Axis          Direction (x y z)        Perp to Cn-Axis"
         DO i=1,NC2PRP
            WRITE (LUW,'(1X,I4,5X,3F10.5,11X,I4)')
     !i,(C2PRP(i,j),j=1,3),IPRPCN
         END DO
      END IF

ccccccccccccccccccc
c Reflexion Planes
ccccccccccccccccccc
      MAXSIGMA=MAXOPSYM
      IF ( SPHERICAL ) THEN
         NC2=0
         DO I=1,NCNAXIS
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(1.D0/2.D0)) NC2=NC2+1
         END DO
         IF      ( NC2.EQ.3 ) THEN
            MAXSIGMA = 6
         ELSE IF ( NC2.EQ.9 ) THEN
            MAXSIGMA = 9
         ELSE IF ( NC2.EQ.15 ) THEN
            MAXSIGMA = 15
         END IF
      END IF
      CALL FINDPLANES (LUW,NSIGMA,SIGMAVECS,MAXSIGMA,Thres)

cccccccccccccccccccccccccccccc
c Improper Rotation Axes, Sn's
cccccccccccccccccccccccccccccc
      NSN=0
      DO i=1,NCNAXIS
         
c First Sn axes
         IF ( (CNAXIS(i,5)/CNAXIS(i,4)).NE.0.5D0 ) THEN
            THETA = 2.D0*PI*CNAXIS(i,5)/CNAXIS(i,4)
            vec(1)=CNAXIS(i,1)
            vec(2)=CNAXIS(i,2)
            vec(3)=CNAXIS(i,3)
            CALL SNMATRIX (vec,SMTZ,THETA)
            CALL CHECKSYMOPER (SMTZ,Thres,CONFIRM)
            IF ( CONFIRM ) 
     !CALL NEWCN(vec,CNAXIS(i,4),CNAXIS(i,5),NSN,SNAXIS)
         END IF

c Now S2n axes
         NMAX=NINT(CNAXIS(i,4))*2
         DO ISN=1,NMAX-1
            THETA = 2.D0*PI*DFLOAT(ISN)/DFLOAT(NMAX)
            vec(1)=CNAXIS(i,1)
            vec(2)=CNAXIS(i,2)
            vec(3)=CNAXIS(i,3)
            CONFIRM=.FALSE.
            IF ( DFLOAT(ISN)/DFLOAT(NMAX).NE.0.5D0 ) THEN
               CALL SNMATRIX (vec,SMTZ,THETA)
               CALL CHECKSYMOPER (SMTZ,Thres,CONFIRM)
            END IF
            IF ( CONFIRM ) 
     !CALL NEWCN(vec,DFLOAT(NMAX),DFLOAT(ISN),NSN,SNAXIS)
         END DO

      END DO
c Printing Sn's found
      IF ( NSN.GT.0 ) THEN
         WRITE (LUW,'(/A,I3/A/A)')
     !" * The number of S_n^k axes found is :",NSN,
     !"   --------------------------------------",
     !" Axis   n   k   theta         Direction (x y z)"
         DO i=1,NSN
            N = NINT(SNAXIS(i,4))
            K = NINT(SNAXIS(i,5))
            THETA = 360.D0*SNAXIS(i,5)/SNAXIS(i,4)
            WRITE (LUW,'(I4,1X,2I4,4X,I4,1X,3F10.5)')
     !i,N,K,NINT(theta),(SNAXIS(i,j),j=1,3)
         END DO
      END IF

c Finding the point group
      IF ( NCNAXIS.GT.0 ) THEN
         ICNMAX=0
         cnx=CNAXIS(1,1)
         cny=CNAXIS(1,2)
         cnz=CNAXIS(1,3)
         DO I=1,NCNAXIS
            ICURR=NINT(CNAXIS(I,4))
            IF ( ICNMAX.LT.ICURR ) THEN
               ICNMAX=ICURR
               cnx=CNAXIS(I,1)
               cny=CNAXIS(I,2)
               cnz=CNAXIS(I,3)
            END IF
         END DO
         iC2PP = 0
         DO I=1,NCNAXIS
            ICURR=NINT(CNAXIS(I,4))
            IF ( ICURR.EQ.2 ) THEN
               dotprod = cnx*CNAXIS(I,1) + cny*CNAXIS(I,2) +
     *                   cnz*CNAXIS(I,3)
               IF ( DABS(dotprod).LT.Thres ) iC2PP = 1
            END IF
         END DO
         ISNMAX=0
         DO I=1,NSN
            ICURR=NINT(SNAXIS(I,4))
            IF ( ISNMAX.LT.ICURR ) ISNMAX=ICURR
         END DO
         SIGMAV=.FALSE.
         IF ( NSIGMA.GT.0 ) THEN
            DO I=1,NSIGMA
               DOT = cnx*SIGMAVECS(i,1) + cny*SIGMAVECS(i,2) + 
     !               cnz*SIGMAVECS(i,3)
               IF ( DABS(DOT).LT.Thres ) SIGMAV=.TRUE.
            END DO
         END IF

         IF ( NC2PRP.GT.0 .OR. iC2PP.EQ.1 ) THEN
c S2n
            IF ( ISNMAX.EQ.2*ICNMAX ) THEN
               IF ( ICNMAX.LT.10 ) WRITE (PG,"('D',I1,'d ')") ICNMAX
               IF ( ICNMAX.GT.10 ) WRITE (PG,"('D',I2,'d')") ICNMAX
c Sn
            ELSE IF ( ISNMAX.EQ.ICNMAX ) THEN
               IF ( ICNMAX.LT.10 ) WRITE (PG,"('D',I1,'h ')") ICNMAX
               IF ( ICNMAX.GT.10 ) WRITE (PG,"('D',I2,'h')") ICNMAX
c INVERSION = S2
            ELSE IF ( INVERSION ) THEN
               IF ( ICNMAX.LT.10 ) WRITE (PG,"('D',I1,'h ')") ICNMAX
               IF ( ICNMAX.GT.10 ) WRITE (PG,"('D',I2,'h')") ICNMAX
            ELSE
               IF ( ICNMAX.LT.10 ) WRITE (PG,"('D',I1,'  ')") ICNMAX
               IF ( ICNMAX.GT.10 ) WRITE (PG,"('D',I2,'  ')") ICNMAX
            END IF
         ELSE
c S2n
            IF ( ISNMAX.EQ.2*ICNMAX ) THEN
               IF ( ICNMAX.LT.10 ) WRITE (PG,"('S',I1,' ')") ISNMAX
               IF ( ICNMAX.GT.10 ) WRITE (PG,"('D',I2,'n')") ICNMAX
c Sn
            ELSE IF ( ISNMAX.EQ.ICNMAX ) THEN
               IF ( ICNMAX.LT.10 ) WRITE (PG,"('C',I1,'h ')") ICNMAX
               IF ( ICNMAX.GT.10 ) WRITE (PG,"('C',I2,'h')") ICNMAX
c INVERSION = S2
            ELSE IF ( INVERSION ) THEN
               IF ( ICNMAX.LT.10 ) WRITE (PG,"('C',I1,'h ')") ICNMAX
               IF ( ICNMAX.GT.10 ) WRITE (PG,"('C',I2,'h')") ICNMAX
c Sigma_v
            ELSE IF ( SIGMAV ) THEN
               IF ( ICNMAX.LT.10 ) WRITE (PG,"('C',I1,'v ')") ICNMAX
               IF ( ICNMAX.GT.10 ) WRITE (PG,"('C',I2,'v')") ICNMAX
            ELSE 
               IF ( ICNMAX.LT.10 ) WRITE (PG,"('C',I1,'  ')") ICNMAX
               IF ( ICNMAX.GT.10 ) WRITE (PG,"('C',I2,' ')") ICNMAX
            END IF
         END IF

      ELSE

         PG = "C1  "
         IF ( NSIGMA.GT.0 ) PG = "Cs  "
         IF ( INVERSION   ) PG = "Ci  "

      END IF

      IF ( SPHERICAL ) THEN
c Getting the number of C2, C3, C4 and C5
         NC2=0
         NC3=0
         NC4=0
         NC5=0
         DO I=1,NCNAXIS
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(1.D0/2.D0)) NC2=NC2+1
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(1.D0/3.D0)) NC3=NC3+1
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(2.D0/3.D0)) NC3=NC3+1
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(1.D0/4.D0)) NC4=NC4+1
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(3.D0/4.D0)) NC4=NC4+1
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(1.D0/5.D0)) NC5=NC5+1
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(2.D0/5.D0)) NC5=NC5+1
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(3.D0/5.D0)) NC5=NC5+1
            IF ( CNAXIS(I,5)/CNAXIS(I,4).EQ.(4.D0/5.D0)) NC5=NC5+1
         END DO
         IF ( NC2.EQ.3 .AND. INVERSION                     ) PG = "Th "
         IF ( NC2.EQ.3 .AND..NOT.INVERSION.AND.NSIGMA.GT.0 ) PG = "Td "
         IF ( NC2.EQ.3 .AND..NOT.INVERSION.AND.NSIGMA.EQ.0 ) PG = "T  "
         IF ( NC2.EQ.9 .AND. INVERSION                     ) PG = "Oh "
         IF ( NC2.EQ.9 .AND..NOT.INVERSION.AND.NSIGMA.GT.0 ) PG = "Od "
         IF ( NC2.EQ.9 .AND..NOT.INVERSION.AND.NSIGMA.EQ.0 ) PG = "O  "
         IF ( NC2.EQ.15.AND. INVERSION                     ) PG = "Ih "
         IF ( NC2.EQ.15.AND..NOT.INVERSION.AND.NSIGMA.GT.0 ) PG = "Id "
         IF ( NC2.EQ.15.AND..NOT.INVERSION.AND.NSIGMA.EQ.0 ) PG = "I  "
      END IF
      CALL WRTPG (LUW,PG)
      NORDPG=NSIGMA+NSN+NCNAXIS+NC2PRP+1
      IF ( INVERSION ) NORDPG=NORDPG+1
      WRITE (LUW,'(A,I4//A,I4)')
     !" * Order of the group: ",NORDPG,
     !" * The number of atoms is: ",NATOM

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine to find reflexion planes
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FINDPLANES(LUW,NSIGMA,SIGMAVECS,MAXSIGMA,Thres)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DIMENSION SIGMAVECS(MAXOPSYM,3), SMTZ(3,3), vec(3)
      LOGICAL CONFIRM

      PI      = 3.14159265358979323846D0

      NSIGMA=0
      DO i=1,nseaset
         IF ( NSEA(i).GT.1 ) THEN
            DO iA = 1,(NSEA(i)-1)
               iatmA = ISEA(i,iA)
               DO iB = (iA+1),NSEA(i)
                  iatmB = ISEA(i,iB)
                  DO j=1,3
                     vec(j) = XYZ(iatmA,j) - XYZ(iatmB,j)
                  END DO
                  vecnorm = DiNORM(vec,3,3)
                  DO j=1,3
                     vec(j) = vec(j) / vecnorm
                  END DO
                  CONFIRM = .FALSE.
                  IF ( NSIGMA.LT.MAXSIGMA ) THEN
                     CALL SIGMA (vec,SMTZ)
                     CALL CHECKSYMOPER (SMTZ,Thres,CONFIRM)
                  END IF
                  IF ( CONFIRM ) THEN
                     IF ( NSIGMA.EQ.0 ) THEN
                        NSIGMA=1
                        DO j=1,3
                           SIGMAVECS(1,j) = vec(j)
                        END DO
                     ELSE
                        IPAR=0
                        DO isigma=1,NSIGMA
                           DOT = vec(1)*SIGMAVECS(isigma,1) +
     *                           vec(2)*SIGMAVECS(isigma,2) +
     *                           vec(3)*SIGMAVECS(isigma,3) 
c Two planes are considered paralell if they normal vectors have an
c angle of 1 degree or less
                           IF ( DABS(DOT).GE.DCOS(PI/180.D0) ) IPAR=1
                        END DO
                        IF ( IPAR.EQ.0 ) THEN
                           NSIGMA=NSIGMA+1
                           IF ( NSIGMA.GT.MAXOPSYM ) THEN
                              WRITE (*,'(//A,I4/A/A,I4//)')
     !" Error: The number of symmetry planes found : ", NSIGMA,
     ! "       is greater than the maximum order allowed for a ",
     ! "       point group ", MAXOPSYM
                              STOP " PLACZEK will stops here "
                           END IF
                           DO j=1,3
                              SIGMAVECS(NSIGMA,j) = vec(j) 
                           END DO
                        END IF
                     END IF
                  END IF
               END DO
            END DO
         END IF
      END DO
c Molecule plane (if any)
      CALL MOLPLANE (NSIGMA,SIGMAVECS)
c printing normal vectors found
      IF ( NSIGMA.GT.0 ) THEN
         WRITE (LUW,'(/A,I3/A/A)')
     !" * The number of reflexion planes found is :",NSIGMA,
     !"   --------------------------------------------",
     !" Plane    Normal vector (x y z)"
         DO i=1,NSIGMA
            WRITE (LUW,'(I4,3X,3F10.5)')i,(SIGMAVECS(i,j),j=1,3)
         END DO
      END IF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C RETURN ATOMIC SYMBOL AND CHARGE
C Input: ATOMIC SYMBOL
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      CHARACTER*2 FUNCTION ATOMNAME( SYMBOL, Z)
      IMPLICIT REAL*8 (A-H,O-Z)

      PARAMETER ( LASTATOM = 94 )
      CHARACTER AN(LASTATOM)*2, SYMBOL*(*)

      ZATM = 0.D0

      DATA (AN(I),I=1,LASTATOM) /
     !"H ","He",
     !"Li","Be","B ","C ","N ","O ","F ","Ne",
     !"Na","Mg","Al","Si","P ","S ","Cl","Ar",
     !"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn",
     !    "Ga","Ge","As","Se", "Br","Kr",
     !"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
     !     "In","Sn","Sb","Te","I ","Xe",
     !"Cs","Ba",
     !     "La","Ce","Pr","Nd","Pm","Sm","Eu",
     !                    "Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
     !          "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
     !     "Tl","Pb","Bi","Po","At","Rn",
     !"Fr","Ra","Ac","Th","Pa","U ","Np","Pu" /

      ATOMNAME="xx"
      DO I=1,LASTATOM
         IF (SYMBOL.EQ.AN(I) ) ATOMNAME = AN(I)
         IF (SYMBOL.EQ.AN(I) ) Z = DFLOAT(I)
      END DO
      IF ( ATOMNAME.EQ."xx" ) THEN
         WRITE (*,"(//A,A)")
     *" Unknown symbol found in input file :",SYMBOL,
     *" Valid symbols are:"
         WRITE (*,"(20(1X,A2))") (AN(I),I= 1,20)
         WRITE (*,"(20(1X,A2))") (AN(I),I=21,40)
         WRITE (*,"(20(1X,A2))") (AN(I),I=41,60)
         WRITE (*,"(20(1X,A2))") (AN(I),I=61,80)
         WRITE (*,"(20(1X,A2))") (AN(I),I=81,94)
         WRITE (*,*)
         STOP
      END IF
  
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION DiNORM(Di,NATOM,MAXAT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Di(MAXAT)

      DiNORM = 0.D0
      DO i=1,NATOM
         DiNORM = DiNORM + Di(i)*Di(i)
      END DO
      DiNORM = DSQRT(DiNORM)    

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Unit matrix of a reflexion plane
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SIGMA(vec,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION vec(3),R(3,3)

      x=vec(1)
      y=vec(2)
      z=vec(3)

      R(1,1) = 1.D0-2.D0*x*x
      R(2,2) = 1.D0-2.D0*y*y
      R(3,3) = 1.D0-2.D0*z*z
      R(1,2) =     -2.D0*x*y
      R(2,1) = R(1,2)
      R(1,3) =     -2.D0*x*z
      R(3,1) = R(1,3)
      R(2,3) =     -2.D0*y*z;
      R(3,2) = R(2,3)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Unit matrix corresponding to a proper rotation
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CNMATRIX(vec,R,THETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION vec(3),R(3,3)

      x=vec(1)
      y=vec(2)
      z=vec(3)

      C=DCOS(THETA)
      S=DSIN(THETA)
      OC=1.D0-C

      R(1,1) = x*x*OC+C
      R(2,2) = y*y*OC+C
      R(3,3) = z*z*OC+C
      R(1,2) = x*y*OC-z*S
      R(2,1) = x*y*OC+z*S
      R(1,3) = x*z*OC+y*S
      R(3,1) = x*z*OC-y*S
      R(2,3) = y*z*OC-x*S
      R(3,2) = y*z*OC+x*S

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Unit matrix corresponding to an improper rotation
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SNMATRIX(vec,SN,THETA)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION vec(3),R(3,3),SIG(3,3),SN(3,3)

      x=vec(1)
      y=vec(2)
      z=vec(3)

      C=DCOS(THETA)
      S=DSIN(THETA)
      OC=1.D0-C

      R(1,1) = x*x*OC+C
      R(2,2) = y*y*OC+C
      R(3,3) = z*z*OC+C
      R(1,2) = x*y*OC-z*S
      R(2,1) = x*y*OC+z*S
      R(1,3) = x*z*OC+y*S
      R(3,1) = x*z*OC-y*S
      R(2,3) = y*z*OC-x*S
      R(3,2) = y*z*OC+x*S

      SIG(1,1) = 1.D0-2.D0*x*x
      SIG(2,2) = 1.D0-2.D0*y*y
      SIG(3,3) = 1.D0-2.D0*z*z
      SIG(1,2) =     -2.D0*x*y
      SIG(2,1) = SIG(1,2)
      SIG(1,3) =     -2.D0*x*z
      SIG(3,1) = SIG(1,3)
      SIG(2,3) =     -2.D0*y*z;
      SIG(3,2) = SIG(2,3)

      CALL PMATRIX (SIG,R,SN,3,3,3)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Verify if that symmetry operation represented by matrix R 
C indeed exists in the molecule
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CHECKSYMOPER(R,Thres,CONFIRM)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION R(3,3), RX(3,MAXAT)
      LOGICAL JUMP(MAXAT), CONFIRM

      CONFIRM = .FALSE.

      DO j=1,NATOM
         DO i=1,3
            RX(i,j) = 0.d0
            DO k=1,3
               RX(i,j) = RX(i,j) + R(i,k) * XYZ(j,k) 
            END DO
         END DO
      END DO

      NEQATM=0
      DO i=1,NSEASET
         DO j=1,NATOM
            JUMP(j) = .FALSE.
         END DO
         DO j=1,NSEA(i)
            iSU = ISEA(i,j)
            XSU = XYZ(iSU,1)
            YSU = XYZ(iSU,2)
            ZSU = XYZ(iSU,3)
            DO j2=1,NSEA(i)
               IF ( .NOT.JUMP(j2) ) THEN
                  iSE = ISEA(i,j2)
                  XSE = RX(1,iSE)
                  YSE = RX(2,iSE)
                  ZSE = RX(3,iSE)
                  DIST =  (XSU-XSE)**2 + (YSU-YSE)**2 + (ZSU-ZSE)**2
                  IF ( DSQRT(DIST).LE.Thres ) THEN
                     NEQATM=NEQATM+1
                     JUMP(j2) = .TRUE.
                  END IF
               END IF
            END DO
         END DO
      END DO
      IF ( NEQATM.EQ.NATOM ) CONFIRM = .TRUE.

      RETURN
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CHECK IF THE MOLECULE IS PLANAR 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MOLPLANE(NSIGMA,SIGMAVECS)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION Tn(3,3), D(3), PMOM(3), SIGMAVECS(MAXOPSYM,3)

      BOHR    = 0.529177208D0
C STARTING EVALUATION OF INERTIA TENSOR

      DO I = 1,3
         DO J = 1,3
           Tn(I,J) = 0.D0
         END DO
      END DO

C EVALUATING CENTER OF MASS COORDINATES: CMX,CMY,CMZ
      CMX   = 0.D0
      CMY   = 0.D0
      CMZ   = 0.D0
      SMASS = 0.D0
      DO I = 1,NATOM
         IF (AU) THEN
            XYZ(I,1) = XYZ(I,1) * BOHR
            XYZ(I,2) = XYZ(I,2) * BOHR
            XYZ(I,3) = XYZ(I,3) * BOHR
         END IF
         AMASS = RMASS(I)
         CMX = CMX + AMASS * XYZ(I,1)
         CMY = CMY + AMASS * XYZ(I,2)
         CMZ = CMZ + AMASS * XYZ(I,3)
         SMASS = SMASS + AMASS
      END DO
      CMX = CMX / SMASS
      CMY = CMY / SMASS
      CMZ = CMZ / SMASS


C INERTIA TENSOR
      DO I = 1,NATOM
         AMASS = RMASS(I)
         CX = XYZ(I,1) - CMX
         CY = XYZ(I,2) - CMY
         CZ = XYZ(I,3) - CMZ

         Tn(1,1) = Tn(1,1) + AMASS * ( CY**2 + CZ**2 )
         Tn(1,2) = Tn(1,2) - AMASS *   CX*CY
         Tn(1,3) = Tn(1,3) - AMASS *   CX*CZ

         Tn(2,1) = Tn(1,2)
         Tn(2,2) = Tn(2,2) + AMASS * ( CX**2 + CZ**2 )
         Tn(2,3) = Tn(2,3) - AMASS *   CY*CZ

         Tn(3,1) = Tn(1,3)
         Tn(3,2) = Tn(2,3)
         Tn(3,3) = Tn(3,3) + AMASS * ( CX**2 + CY**2 )
      END DO

C DIAGONALIZATION OF INERTIA TENSOR
      CALL JACOB (Tn,D,3,3)

      DO I=1,3
         PMOM(I) = D(I)
      END DO
      CALL REORD (Tn,D,3,3)

C SETTING IA,IB,IC
      dIC = D(1)
      dIB = D(2)
      dIA = D(3)

C ROTATIONAL SYMMETRY
      t1 = DABS( dIA - dIB )
      t2 = DABS( dIA - dIC )
      t3 = DABS( dIB - dIC )

C SEARCHING FOR LINEARITY & PLANARITY
      Thres = 1.D-04
      LINEAR=.FALSE.
      IF ( DABS(dIA).LT.Thres .AND. DABS(dIB-dIC).LT.Thres )
     !LINEAR = .TRUE.
      IF ( ABS(dIC-dIB-dIA).LT.Thres .AND. .NOT.LINEAR ) THEN
         NSIGMA=NSIGMA+1
c The normal vector of the molecular plane is IC
         DO i=1,3
            SIGMAVECS(NSIGMA,i) = Tn(i,1)
         END DO
      END IF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS ROUTINE CALCULATES THE PRINCIPAL MOMENTS OF SEA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SEAxIAxIBxIC(iseaset,PIA,PIB,PIC,VA,VB,VC,CMX,CMY,CMZ)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION Tn(3,3), D(3),VA(3),VB(3),VC(3)

c Calculating the Center of Mass of the SEA set
      CMX = 0.D0
      CMY = 0.D0
      CMZ = 0.D0
      TMASS = 0.D0
      DO j = 1,NSEA(iseaset)
         AMASS = RMASS(ISEA(iseaset,j))
         CMX = CMX + XYZ(ISEA(iseaset,j),1) * AMASS
         CMY = CMY + XYZ(ISEA(iseaset,j),2) * AMASS
         CMZ = CMZ + XYZ(ISEA(iseaset,j),3) * AMASS
         TMASS = TMASS + AMASS
      END DO
      CMX = CMX / TMASS
      CMY = CMY / TMASS
      CMZ = CMZ / TMASS

C INERTIA TENSOR
      DO I = 1,3
         DO J = 1,3
           Tn(I,J) = 0.D0
         END DO
      END DO
      DO j = 1,NSEA(iseaset)
         CX = XYZ(ISEA(iseaset,j),1)-CMX
         CY = XYZ(ISEA(iseaset,j),2)-CMY
         CZ = XYZ(ISEA(iseaset,j),3)-CMZ 

         AMASS = RMASS(ISEA(iseaset,j))

         Tn(1,1) = Tn(1,1) + AMASS * ( CY**2 + CZ**2 )
         Tn(1,2) = Tn(1,2) - AMASS *   CX*CY
         Tn(1,3) = Tn(1,3) - AMASS *   CX*CZ
         Tn(2,1) = Tn(1,2)
         Tn(2,2) = Tn(2,2) + AMASS * ( CX**2 + CZ**2 )
         Tn(2,3) = Tn(2,3) - AMASS *   CY*CZ
         Tn(3,1) = Tn(1,3)
         Tn(3,2) = Tn(2,3)
         Tn(3,3) = Tn(3,3) + AMASS * ( CX**2 + CY**2 )
      END DO

C DIAGONALIZATION OF INERTIA TENSOR
      CALL JACOB (Tn,D,3,3)
      CALL REORD (Tn,D,3,3)

C SETTING IA,IB,IC
      PIC = D(1)
      PIB = D(2)
      PIA = D(3)

      DO i=1,3
        VA(i)=Tn(i,3)
        VB(i)=Tn(i,2)
        VC(i)=Tn(i,1)
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS ROUTINE CALCULATES THE PRINCIPAL MOMENTS OF THE ENTIRE MOLECULE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MOLxIAxIBxIC(PIA,PIB,PIC,VA,VB,VC)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION Tn(3,3), D(3),VA(3),VB(3),VC(3)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Calculating the Center of Mass 
      CMX = 0.D0
      CMY = 0.D0
      CMZ = 0.D0
      TMASS = 0.D0
      DO j = 1,NATOM
         AMASS = RMASS(j)
         CMX = CMX + XYZ(j,1) * AMASS
         CMY = CMY + XYZ(j,2) * AMASS
         CMZ = CMZ + XYZ(j,3) * AMASS
         TMASS = TMASS + AMASS
      END DO
      CMX = CMX / TMASS
      CMY = CMY / TMASS
      CMZ = CMZ / TMASS

C INERTIA TENSOR
      DO I = 1,3
         DO J = 1,3
           Tn(I,J) = 0.D0
         END DO
      END DO
      DO j = 1,NATOM
         CX = XYZ(j,1)-CMX
         CY = XYZ(j,2)-CMY
         CZ = XYZ(j,3)-CMZ 

         AMASS = RMASS(j)

         Tn(1,1) = Tn(1,1) + AMASS * ( CY**2 + CZ**2 )
         Tn(1,2) = Tn(1,2) - AMASS *   CX*CY
         Tn(1,3) = Tn(1,3) - AMASS *   CX*CZ
         Tn(2,1) = Tn(1,2)
         Tn(2,2) = Tn(2,2) + AMASS * ( CX**2 + CZ**2 )
         Tn(2,3) = Tn(2,3) - AMASS *   CY*CZ
         Tn(3,1) = Tn(1,3)
         Tn(3,2) = Tn(2,3)
         Tn(3,3) = Tn(3,3) + AMASS * ( CX**2 + CY**2 )
      END DO

C DIAGONALIZATION OF INERTIA TENSOR
      CALL JACOB (Tn,D,3,3)
      CALL REORD (Tn,D,3,3)

C SETTING IA,IB,IC
      PIC = D(1)
      PIB = D(2)
      PIA = D(3)

      DO i=1,3
        VA(i)=Tn(i,3)
        VB(i)=Tn(i,2)
        VC(i)=Tn(i,1)
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS ROUTINE LOCATE PERPENDICULAR C2 AXES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE C2PERP(NC2PRP,C2PRP,NCNAXIS,CNAXIS,Thres,NSYMVEC,
     *                  SYMVEC,IPRPCN)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION C2PRP(MAXOPSYM,5), CNMTZ(3,3), vec(3), 
     *          CNAXIS(MAXOPSYM,5), 
     *          SYMVEC(MAXOPSYM,3), Va(3), Vb(3)
      LOGICAL CONFIRM

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0

c Identifying the highest order Cn in the molecule
      IPRPCN=0
      ICNMAX=0
      DO K=1,NCNAXIS
         IF ( NINT(CNAXIS(K,4)).GT.ICNMAX ) THEN
            ICNMAX = NINT(CNAXIS(K,4))
            IPRPCN=K
         END IF
      END DO

c Try to find C2 from middle point vectors of two atoms SEA sets
      DO I=1,NSEASET
         IF ( NSEA(I).GT.1 ) THEN

            DO J =1,NSEA(I)-1
               DO J2=J+1,NSEA(I)

      vec(1) = ( XYZ(ISEA(I,J),1)+XYZ(ISEA(I,J2),1) ) / 2.D0
      vec(2) = ( XYZ(ISEA(I,J),2)+XYZ(ISEA(I,J2),2) ) / 2.D0
      vec(3) = ( XYZ(ISEA(I,J),3)+XYZ(ISEA(I,J2),3) ) / 2.D0

                  VECNORM = DSQRT( vec(1)**2 + vec(2)**2 + vec(3)**2 )
                  IF ( VECNORM.GT.Thres ) THEN
                     vec(1) = vec(1) / VECNORM
                     vec(2) = vec(2) / VECNORM
                     vec(3) = vec(3) / VECNORM
                     DOT = vec(1)*CNAXIS(IPRPCN,1) + 
     *                     vec(2)*CNAXIS(IPRPCN,2) + 
     *                     vec(3)*CNAXIS(IPRPCN,3)
                     IF ( DABS(DOT).LT.Thres ) THEN
                        THETA = PI
                        CALL CNMATRIX (vec,CNMTZ,THETA)
                        CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
                        IF ( CONFIRM ) THEN
                           NCNAXISp = NCNAXIS
                           CALL NEWCN (vec,2.D0,1.D0,NCNAXIS,CNAXIS)
                           IF ( (NCNAXIS-1).EQ.NCNAXISp ) THEN
                              CALL NEWCN (vec,2.D0,1.D0,NC2PRP,C2PRP)
                              NCNAXIS = NCNAXISp
                           END IF
c                          IF ( NC2PRP.EQ.NSEA(I) ) RETURN
                        END IF
                     END IF
                  END IF

               END DO
            END DO

         END IF
      END DO

c Try to find C2 from position vector of one atom from the SEA set
      DO I=1,NSEASET
         DO J =1,NSEA(I)

            vec(1) = XYZ(ISEA(I,J),1)
            vec(2) = XYZ(ISEA(I,J),2)
            vec(3) = XYZ(ISEA(I,J),3)

            VECNORM = DSQRT( vec(1)**2 + vec(2)**2 + vec(3)**2 )
            IF ( VECNORM.GT.Thres ) THEN
               vec(1) = vec(1) / VECNORM
               vec(2) = vec(2) / VECNORM
               vec(3) = vec(3) / VECNORM
               DOT = vec(1)*CNAXIS(IPRPCN,1) + 
     *               vec(2)*CNAXIS(IPRPCN,2) + 
     *               vec(3)*CNAXIS(IPRPCN,3)
               IF ( DABS(DOT).LT.Thres ) THEN
                  THETA = PI
                  CALL CNMATRIX (vec,CNMTZ,THETA)
                  CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
                  IF ( CONFIRM ) THEN
                     NCNAXISp = NCNAXIS
                     CALL NEWCN (vec,2.D0,1.D0,NCNAXIS,CNAXIS)
                     IF ( (NCNAXIS-1).EQ.NCNAXISp ) THEN
                        CALL NEWCN (vec,2.D0,1.D0,NC2PRP,C2PRP)
                        NCNAXIS = NCNAXISp
                     END IF
c                    IF ( NC2PRP.EQ.NSEA(I) ) RETURN
                  END IF
               END IF
            END IF

         END DO
      END DO

c Try to find C2 from those vectors connecting two atoms of some SEA set
c using a crossproduct
      IF (NSYMVEC.LE.1) RETURN
      DO J =1,NSYMVEC-1
      DO J2=J,NSYMVEC

         DO K=1,3
            Va(K)=SYMVEC(J,K)
            Vb(K)=SYMVEC(J2,K)
            CALL CROSSPROD (Va,Vb,vec)
         END DO

         VECNORM = DSQRT( vec(1)**2 + vec(2)**2 + vec(3)**2 )
         IF ( VECNORM.GT.Thres ) THEN
            vec(1) = vec(1) / VECNORM
            vec(2) = vec(2) / VECNORM
            vec(3) = vec(3) / VECNORM
            THETA = PI
            CALL CNMATRIX (vec,CNMTZ,THETA)
            CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
            IF ( CONFIRM ) THEN
               DOT = vec(1)*CNAXIS(IPRPCN,1) +
     *               vec(2)*CNAXIS(IPRPCN,2) + 
     *               vec(3)*CNAXIS(IPRPCN,3)
               IF ( DABS(DOT).LT.Thres ) THEN
                  NCNAXISp = NCNAXIS
                  CALL NEWCN (vec,2.D0,1.D0,NCNAXIS,CNAXIS)
                  IF ( (NCNAXIS-1).EQ.NCNAXISp ) THEN
                     CALL NEWCN (vec,2.D0,1.D0,NC2PRP,C2PRP)
                     NCNAXIS = NCNAXISp
                  END IF
               END IF
            END IF
         END IF

      END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS ROUTINE WILL CONFIRM IF A GIVEN CN EXIST IN THE MOLECULE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CNCONFIRM(NMAX,VCN,Thres,NCNAXIS,CNAXIS)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION CNMTZ(3,3), VCN(3), CNAXIS(MAXOPSYM ,5)
      LOGICAL CONFIRM

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0

      DO ICN=1,NMAX
         IF ( ICN/NMAX.GE.1 ) RETURN
         THETA = 2.D0*PI*DFLOAT(ICN)/DFLOAT(NMAX)
         CALL CNMATRIX (VCN,CNMTZ,THETA)
         CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
         IF ( CONFIRM ) THEN
c Este procedimento descarta o novo C2 para ele ser encontrado pela
c rotina do C2 perpendicular.
c            NOTCONF=0
c            DO K=1,NCNAXIS
c               DOT = VCN(1)*CNAXIS(K,1) +
c     *               VCN(2)*CNAXIS(K,2) + 
c     *               VCN(3)*CNAXIS(K,3)
cc                  THETAp=2.D0*PI*CNAXIS(K,5)/CNAXIS(K,4)
c                  IF ( DABS(DOT).LT.Thres .AND.
c     *                 DABS(THETA-PI).LT.Thres ) NOTCONF=1
c            END DO
c            IF (NOTCONF.EQ.0)
            CALL NEWCN (VCN,DFLOAT(NMAX),DFLOAT(ICN),NCNAXIS,CNAXIS)
         END IF
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS ROUTINE LOCATE C2 AXES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE C2AXIS(NSYMVEC,SYMVEC,I,NCNAXIS,CNAXIS,Thres)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION CNAXIS(MAXOPSYM,5), CNMTZ(3,3), V1(3),V2(3),vec(3),
     !          SYMVEC(MAXOPSYM,3)

      LOGICAL CONFIRM

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0

c C2 from middle points
      DO J =1,NSEA(I)-1
         DO J2 =J+1,NSEA(I)
            vec(1) = ( XYZ(ISEA(I,J),1)+XYZ(ISEA(I,J2),1) ) / 2.D0
            vec(2) = ( XYZ(ISEA(I,J),2)+XYZ(ISEA(I,J2),2) ) / 2.D0
            vec(3) = ( XYZ(ISEA(I,J),3)+XYZ(ISEA(I,J2),3) ) / 2.D0

            VECNORM = DSQRT( vec(1)**2 + vec(2)**2 + vec(3)**2 )
            IF ( VECNORM.GT.Thres ) THEN
               vec(1) = vec(1) / VECNORM
               vec(2) = vec(2) / VECNORM
               vec(3) = vec(3) / VECNORM
               THETA = PI
               CALL CNMATRIX (vec,CNMTZ,THETA)
               CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
               IF ( CONFIRM ) THEN
                  NOTCONF=0
                  DO K=1,NCNAXIS
                     DOT = vec(1)*CNAXIS(K,1) + 
     !                     vec(2)*CNAXIS(K,2) + 
     !                     vec(3)*CNAXIS(K,3)
                     THETAp=2.D0*PI*CNAXIS(K,5)/CNAXIS(K,4)
                     THETA =     PI
                     IF ( DABS(DOT).LT.Thres .AND. 
     !                    DABS(THETA-THETAp).LT.(PI/180.D0) ) NOTCONF=1
                  END DO
                  IF (NOTCONF.EQ.0) 
     !CALL NEWCN(vec,2.D0,1.D0,NCNAXIS,CNAXIS)
               END IF
            END IF
         END DO
      END DO

c Try to find C2 from position vector of one atom from the SEA set
      DO J =1,NSEA(I)

         vec(1) = XYZ(ISEA(I,J),1)
         vec(2) = XYZ(ISEA(I,J),2)
         vec(3) = XYZ(ISEA(I,J),3)

         VECNORM = DSQRT( vec(1)**2 + vec(2)**2 + vec(3)**2 )
         IF ( VECNORM.GT.Thres ) THEN
            vec(1) = vec(1) / VECNORM
            vec(2) = vec(2) / VECNORM
            vec(3) = vec(3) / VECNORM
            THETA = PI
            CALL CNMATRIX (vec,CNMTZ,THETA)
            CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
            IF ( CONFIRM ) THEN
               NOTCONF=0
               DO K=1,NCNAXIS
                  DOT = vec(1)*CNAXIS(K,1) +
     !                  vec(2)*CNAXIS(K,2) + 
     !                  vec(3)*CNAXIS(K,3)
                  THETAp=2.D0*PI*CNAXIS(K,5)/CNAXIS(K,4)
                  THETA =     PI
                  IF ( DABS(DOT).LT.Thres .AND.
     !                 DABS(THETA-THETAp).LT.(PI/180.D0) ) NOTCONF=1
               END DO
               IF (NOTCONF.EQ.0)
     !CALL NEWCN(vec,2.D0,1.D0,NCNAXIS,CNAXIS)
            END IF
         END IF
      END DO

c Trying to find C2 from cross produts of SYMVECs from SEA sets of
c two atoms
      IF (NSYMVEC.EQ.1) RETURN
      DO J =1,NSYMVEC-1
      DO J2=J,NSYMVEC

         DO K=1,3
            V1(K)=SYMVEC(J,K)
            V2(K)=SYMVEC(J2,K)
            CALL CROSSPROD (V1,V2,Vec)
         END DO

         VECNORM = DSQRT( vec(1)**2 + vec(2)**2 + vec(3)**2 )
         IF ( VECNORM.GT.Thres ) THEN
            vec(1) = vec(1) / VECNORM
            vec(2) = vec(2) / VECNORM
            vec(3) = vec(3) / VECNORM
            THETA = PI
            CALL CNMATRIX (vec,CNMTZ,THETA)
            CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
            IF ( CONFIRM ) THEN
               NOTCONF=0
               DO K=1,NCNAXIS
                  DOT = vec(1)*CNAXIS(K,1) + 
     !                  vec(2)*CNAXIS(K,2) + 
     !                  vec(3)*CNAXIS(K,3)
                  THETAp=2.D0*PI*CNAXIS(K,5)/CNAXIS(K,4)
                  IF ( DABS(DOT).LT.Thres .AND.
     !                 DABS(PI-THETAp).LT.(PI/180.D0) ) NOTCONF=1
               END DO
               IF (NOTCONF.EQ.0) 
     !CALL NEWCN(vec,2.D0,1.D0,NCNAXIS,CNAXIS)
            END IF
         END IF

      END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS ROUTINE LOCATE C2 IN SPHERICAL MOLECULES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE C2SPHER(I,NCNAXIS,CNAXIS,Thres,MAXC2,NC2)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION CNAXIS(MAXOPSYM,5), CNMTZ(3,3), vec(3)
      LOGICAL CONFIRM

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0

      IF ( NSEA(I).LT.2 ) RETURN
      IF ( NC2.GE.MAXC2 ) RETURN

c C2 from middle points
      DO J =1,NSEA(I)-1
         DO J2 =J+1,NSEA(I)
            IF ( NC2.GE.MAXC2 ) RETURN
            vec(1) = ( XYZ(ISEA(I,J),1)+XYZ(ISEA(I,J2),1) ) / 2.D0
            vec(2) = ( XYZ(ISEA(I,J),2)+XYZ(ISEA(I,J2),2) ) / 2.D0
            vec(3) = ( XYZ(ISEA(I,J),3)+XYZ(ISEA(I,J2),3) ) / 2.D0

            VECNORM = DSQRT( vec(1)**2 + vec(2)**2 + vec(3)**2 )
            IF ( VECNORM.GT.Thres ) THEN
               vec(1) = vec(1) / VECNORM
               vec(2) = vec(2) / VECNORM
               vec(3) = vec(3) / VECNORM
               THETA = PI
               CALL CNMATRIX (vec,CNMTZ,THETA)
               CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
               IF ( CONFIRM ) THEN
                  NCNAXISp=NCNAXIS
                  CALL NEWCN (vec,2.D0,1.D0,NCNAXIS,CNAXIS)
                  IF (NCNAXIS.GT.NCNAXISp ) NC2=NC2+1
               END IF
            END IF
         END DO
      END DO

c Try to find C2 from position vector of one atom from the SEA set
      DO J =1,NSEA(I)

         IF ( NC2.GE.MAXC2 ) RETURN
         vec(1) = XYZ(ISEA(I,J),1)
         vec(2) = XYZ(ISEA(I,J),2)
         vec(3) = XYZ(ISEA(I,J),3)

         VECNORM = DSQRT( vec(1)**2 + vec(2)**2 + vec(3)**2 )
         IF ( VECNORM.GT.Thres ) THEN
            vec(1) = vec(1) / VECNORM
            vec(2) = vec(2) / VECNORM
            vec(3) = vec(3) / VECNORM
            THETA = PI
            CALL CNMATRIX (vec,CNMTZ,THETA)
            CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
            IF ( CONFIRM ) THEN
               NCNAXISp=NCNAXIS
               CALL NEWCN (vec,2.D0,1.D0,NCNAXIS,CNAXIS)
               IF (NCNAXIS.GT.NCNAXISp ) NC2=NC2+1
            END IF
         END IF
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS ROUTINE LOCATE C3 AXES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE C3AXIS(I,NCNAXIS,CNAXIS,Thres,MAXC3,NC3)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION CNAXIS(MAXOPSYM,5), CNMTZ(3,3), 
     *          vecA(3),vecB(3),vecC(3)
      LOGICAL CONFIRM

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0

      IF ( NSEA(I).LT.3 ) RETURN
      IF ( NC3.GE.MAXC3 ) RETURN

      DO J1=1,NSEA(I)-2
         DO J2=J1+1,NSEA(I)-1
            DO J3=J2+1,NSEA(I)

               IF ( NC3.GE.MAXC3 ) RETURN
               DO K=1,3
                  vecA(K) = ( XYZ(ISEA(I,J1),K)-XYZ(ISEA(I,J2),K) )
                  vecB(K) = ( XYZ(ISEA(I,J1),K)-XYZ(ISEA(I,J3),K) )
                  vecC(K) = ( XYZ(ISEA(I,J2),K)-XYZ(ISEA(I,J3),K) )
               END DO
               RNORMA=DSQRT(vecA(1)**2+vecA(2)**2+vecA(3)**2)
               RNORMB=DSQRT(vecB(1)**2+vecB(2)**2+vecB(3)**2)
               RNORMC=DSQRT(vecC(1)**2+vecC(2)**2+vecC(3)**2)

               IF ( DABS(RNORMA-RNORMB).LT.Thres .AND. 
     !              DABS(RNORMA-RNORMC).LT.Thres .AND.
     !              DABS(RNORMB-RNORMC).LT.Thres ) THEN

                  CALL CROSSPROD (vecA,vecB,vecC)
                  VECNORM = DSQRT(vecC(1)**2+vecC(2)**2+vecC(3)**2)
                  IF ( VECNORM.GT.Thres ) THEN
                     vecC(1) = vecC(1) / VECNORM
                     vecC(2) = vecC(2) / VECNORM
                     vecC(3) = vecC(3) / VECNORM
c C_3 rotation
                     THETA = 2.D0*PI*1.D0/3.D0
                     CALL CNMATRIX (vecC,CNMTZ,THETA)
                     CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
                     IF ( CONFIRM ) THEN
                        NCNAXISp=NCNAXIS
                        CALL NEWCN (vecC,3.D0,1.D0,NCNAXIS,CNAXIS)
                        IF (NCNAXIS.GT.NCNAXISp ) NC3=NC3+1
c C_3^2 rotation
                        THETA = 2.D0*PI*2.D0/3.D0
                        CALL CNMATRIX (vecC,CNMTZ,THETA)
                        CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
                        IF ( CONFIRM ) THEN
                           NCNAXISp=NCNAXIS
                           CALL NEWCN (vecC,3.D0,2.D0,NCNAXIS,CNAXIS)
                           IF (NCNAXIS.GT.NCNAXISp ) NC3=NC3+1
                        END IF
                     END IF
                  END IF
               END IF

            END DO
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS ROUTINE LOCATE C4 AXES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE C4AXIS(I,NCNAXIS,CNAXIS,Thres,MAXC4,NC4)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION CNAXIS(MAXOPSYM,5), CNMTZ(3,3), 
     *          vecA(3),vecB(3),vecC(3),vecD(3),vecE(3)
      LOGICAL CONFIRM, PLANAR, SQUARE

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0

      IF ( NSEA(I).LT.4 ) RETURN
      IF ( NC4.GE.MAXC4 ) RETURN

      DO J1=1,NSEA(I)-3
         DO J2=J1+1,NSEA(I)-2
            DO J3=J2+1,NSEA(I)-1
               DO J4=J3+1,NSEA(I)

               IF ( NC4.GE.MAXC4 ) RETURN

               PLANAR=.FALSE.
               SQUARE=.FALSE.
               DO K=1,3
                  vecA(K) = ( XYZ(ISEA(I,J1),K)-XYZ(ISEA(I,J2),K) )
                  vecB(K) = ( XYZ(ISEA(I,J1),K)-XYZ(ISEA(I,J3),K) )
                  vecC(K) = ( XYZ(ISEA(I,J1),K)-XYZ(ISEA(I,J4),K) )
               END DO
               CALL CROSSPROD (vecA,vecB,vecD)
               CALL CROSSPROD (vecA,vecC,vecE)
               RNORMD=DSQRT(vecD(1)**2+vecD(2)**2+vecD(3)**2)
               RNORME=DSQRT(vecE(1)**2+vecE(2)**2+vecE(3)**2)
               DO K=1,3
                  IF (RNORMD.GT.Thres) vecD(K)=vecD(K)/RNORMD
                  IF (RNORME.GT.Thres) vecE(K)=vecE(K)/RNORME
               END DO
               DOT=vecD(1)*vecE(1)+vecD(2)*vecE(2)+vecD(3)*vecE(3)
               IF ( DABS(DOT).GE.DCOS(PI/180.D0) ) PLANAR=.TRUE.

               IF ( PLANAR ) THEN
                  RNORMA=DSQRT(vecA(1)**2+vecA(2)**2+vecA(3)**2)
                  RNORMB=DSQRT(vecB(1)**2+vecB(2)**2+vecB(3)**2)
                  RNORMC=DSQRT(vecC(1)**2+vecC(2)**2+vecC(3)**2)
                  DO K=1,3
                     IF (RNORMA.GT.Thres) vecA(K)=vecA(K)/RNORMA
                     IF (RNORMB.GT.Thres) vecB(K)=vecB(K)/RNORMB
                     IF (RNORMC.GT.Thres) vecC(K)=vecC(K)/RNORMC
                  END DO
                  DOT90=vecA(1)*vecB(1)+vecA(2)*vecB(2)+vecA(3)*vecB(3)
                  DOT45=vecA(1)*vecC(1)+vecA(2)*vecC(2)+vecA(3)*vecC(3)
                  IF ( DABS(RNORMA-RNORMB).LT.Thres .AND.
     !                 DABS(DSQRT(2.D0)*RNORMA-RNORMC).LT.Thres .AND.
     !                 DABS(DOT90).LT.Thres .AND.
     !                 DABS(DOT45-DCOS(PI/4.D0)).LT.Thres )SQUARE=.TRUE.
                  DOT90=vecA(1)*vecC(1)+vecA(2)*vecC(2)+vecA(3)*vecC(3)
                  DOT45=vecA(1)*vecB(1)+vecA(2)*vecB(2)+vecA(3)*vecB(3)
                  IF ( DABS(RNORMA-RNORMC).LT.Thres .AND.
     !                 DABS(DSQRT(2.D0)*RNORMA-RNORMB).LT.Thres .AND.
     !                 DABS(DOT90).LT.Thres .AND.
     !                 DABS(DOT45-DCOS(PI/4.D0)).LT.Thres )SQUARE=.TRUE.
                  DOT90=vecB(1)*vecC(1)+vecB(2)*vecC(2)+vecB(3)*vecC(3)
                  DOT45=vecB(1)*vecA(1)+vecB(2)*vecA(2)+vecB(3)*vecA(3)
                  IF ( DABS(RNORMB-RNORMC).LT.Thres .AND.
     !                 DABS(DSQRT(2.D0)*RNORMB-RNORMA).LT.Thres .AND.
     !                 DABS(DOT90).LT.Thres .AND.
     !                 DABS(DOT45-DCOS(PI/4.D0)).LT.Thres )SQUARE=.TRUE.
               END IF

               IF ( SQUARE ) THEN
c C_4 rotation
                  THETA = 2.D0*PI*1.D0/4.D0
                  CALL CNMATRIX (vecD,CNMTZ,THETA)
                  CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
                  IF ( CONFIRM ) THEN
                     NCNAXISp=NCNAXIS
                     CALL NEWCN (vecD,4.D0,1.D0,NCNAXIS,CNAXIS)
                     IF (NCNAXIS.GT.NCNAXISp ) NC4=NC4+1
c C_4^3 rotation
                     THETA = 2.D0*PI*3.D0/4.D0
                     CALL CNMATRIX (vecD,CNMTZ,THETA)
                     CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
                     IF ( CONFIRM ) THEN
                        NCNAXISp=NCNAXIS
                        CALL NEWCN (vecD,4.D0,3.D0,NCNAXIS,CNAXIS)
                        IF (NCNAXIS.GT.NCNAXISp ) NC4=NC4+1
                     END IF
                  END IF
               END IF

               END DO
            END DO
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C THIS ROUTINE LOCATE C5 AXES
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE C5AXIS(I,NCNAXIS,CNAXIS,Thres,MAXC5,NC5)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION CNAXIS(MAXOPSYM,5), CNMTZ(3,3), vec(3),
     *          vecA(3),vecB(3),vecC(3),vecD(3),vecE(3),vecF(3),vecG(3)
      LOGICAL CONFIRM, PLANAR, PENTAGON

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0

      IF ( NSEA(I).LT.5 ) RETURN
      IF ( NC5.GE.MAXC5 ) RETURN

      DO J1=1,NSEA(I)-4
         DO J2=J1+1,NSEA(I)-3
            DO J3=J2+1,NSEA(I)-2
               DO J4=J3+1,NSEA(I)-1
                  DO J5=J4+1,NSEA(I)

      IF ( NC5.GE.MAXC5 ) RETURN

      PLANAR=.FALSE.
      PENTAGON=.FALSE.
      DO K=1,3
         vecA(K) = ( XYZ(ISEA(I,J1),K)-XYZ(ISEA(I,J2),K) )
         vecB(K) = ( XYZ(ISEA(I,J1),K)-XYZ(ISEA(I,J3),K) )
         vecC(K) = ( XYZ(ISEA(I,J1),K)-XYZ(ISEA(I,J4),K) )
         vecD(K) = ( XYZ(ISEA(I,J1),K)-XYZ(ISEA(I,J5),K) )
      END DO
      CALL CROSSPROD (vecA,vecB,vecE)
      CALL CROSSPROD (vecA,vecC,vecF)
      CALL CROSSPROD (vecA,vecD,vecG)
      RNORME=DSQRT(vecE(1)**2+vecE(2)**2+vecE(3)**2)
      RNORMF=DSQRT(vecF(1)**2+vecF(2)**2+vecF(3)**2)
      RNORMG=DSQRT(vecG(1)**2+vecG(2)**2+vecG(3)**2)
      DO K=1,3
         IF (RNORME.GT.Thres) vecE(K)=vecE(K)/RNORME
         IF (RNORMF.GT.Thres) vecF(K)=vecF(K)/RNORMF
         IF (RNORMG.GT.Thres) vecG(K)=vecG(K)/RNORMG
      END DO
      DOT1=vecE(1)*vecF(1)+vecE(2)*vecF(2)+vecE(3)*vecF(3)
      DOT2=vecE(1)*vecG(1)+vecE(2)*vecG(2)+vecE(3)*vecG(3)
      IF ( DABS(DOT1).GE.DCOS(PI/180.D0) .AND.
     !     DABS(DOT2).GE.DCOS(PI/180.D0) ) PLANAR=.TRUE.

      IF ( PLANAR ) THEN
         RNORMA=DSQRT(vecA(1)**2+vecA(2)**2+vecA(3)**2)
         RNORMB=DSQRT(vecB(1)**2+vecB(2)**2+vecB(3)**2)
         RNORMC=DSQRT(vecC(1)**2+vecC(2)**2+vecC(3)**2)
         RNORMD=DSQRT(vecD(1)**2+vecD(2)**2+vecD(3)**2)
         DO K=1,3
            IF (RNORMA.GT.Thres) vecA(K)=vecA(K)/RNORMA
            IF (RNORMB.GT.Thres) vecB(K)=vecB(K)/RNORMB
            IF (RNORMC.GT.Thres) vecC(K)=vecC(K)/RNORMC
            IF (RNORMD.GT.Thres) vecD(K)=vecD(K)/RNORMD
         END DO

         DOTAB=vecA(1)*vecB(1)+vecA(2)*vecB(2)+vecA(3)*vecB(3)
         DOTAC=vecA(1)*vecC(1)+vecA(2)*vecC(2)+vecA(3)*vecC(3)
         DOTAD=vecA(1)*vecD(1)+vecA(2)*vecD(2)+vecA(3)*vecD(3)
         pentagon=.true.
c         IF      ( DABS(DOTAB-DCOS(PI*0.6D0)).LT.Thres ) THEN
c               IF ( DABS(DOTAC-DCOS(0.2D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.4D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c               IF ( DABS(DOTAC-DCOS(0.4D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.2D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c         ELSE IF ( DABS(DOTAB-DCOS(PI*0.4D0)).LT.Thres ) THEN
c               IF ( DABS(DOTAC-DCOS(0.2D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.6D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c               IF ( DABS(DOTAC-DCOS(0.6D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.2D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c         ELSE IF ( DABS(DOTAB-DCOS(PI*0.2D0)).LT.Thres ) THEN
c               IF ( DABS(DOTAC-DCOS(0.4D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.6D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c               IF ( DABS(DOTAC-DCOS(0.6D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.4D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c         ELSE IF ( DABS(DOTAB-DCOS(PI*0.2D0)).LT.Thres ) THEN
c               IF ( DABS(DOTAC-DCOS(0.2D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.4D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c               IF ( DABS(DOTAC-DCOS(0.4D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.2D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c         ELSE IF ( DABS(DOTAB-DCOS(PI*0.4D0)).LT.Thres ) THEN
c               IF ( DABS(DOTAC-DCOS(0.2D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.6D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c               IF ( DABS(DOTAC-DCOS(0.6D0*PI)).LT.Thres .AND.
c     !              DABS(DOTAD-DCOS(0.2D0*PI)).LT.Thres )
c     !            PENTAGON=.TRUE.
c         END IF

         IF ( PENTAGON ) THEN
            DO K=1,3
               IF ( RNORME.GT.Thres ) THEN
                  vec(K) = vecE(K)
               ELSE IF ( RNORMF.GT.Thres ) THEN
                  vec(K) = vecF(K)
               ELSE
                  vec(K) = vecG(K)
               END IF
            END DO
c C_5 rotation
            THETA = 2.D0*PI*1.D0/5.D0
            CALL CNMATRIX (vec,CNMTZ,THETA)
            CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
            IF ( CONFIRM ) THEN
                 NCNAXISp=NCNAXIS
                 CALL NEWCN (vec,5.D0,1.D0,NCNAXIS,CNAXIS)
                 IF (NCNAXIS.GT.NCNAXISp ) NC5=NC5+1
            END IF
c C_5^2 rotation
            THETA = 2.D0*PI*2.D0/5.D0
            CALL CNMATRIX (vec,CNMTZ,THETA)
            CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
            IF ( CONFIRM ) THEN
               NCNAXISp=NCNAXIS
               CALL NEWCN (vec,5.D0,2.D0,NCNAXIS,CNAXIS)
               IF (NCNAXIS.GT.NCNAXISp ) NC5=NC5+1
            END IF
c C_5^3 rotation
            THETA = 2.D0*PI*3.D0/5.D0
            CALL CNMATRIX (vec,CNMTZ,THETA)
            CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
            IF ( CONFIRM ) THEN
               NCNAXISp=NCNAXIS
               CALL NEWCN (vec,5.D0,3.D0,NCNAXIS,CNAXIS)
               IF (NCNAXIS.GT.NCNAXISp ) NC5=NC5+1
            END IF
c C_5^4 rotation
            THETA = 2.D0*PI*4.D0/5.D0
            CALL CNMATRIX (vec,CNMTZ,THETA)
            CALL CHECKSYMOPER (CNMTZ,Thres,CONFIRM)
            IF ( CONFIRM ) THEN
               NCNAXISp=NCNAXIS
               CALL NEWCN (vec,5.D0,4.D0,NCNAXIS,CNAXIS)
               IF (NCNAXIS.GT.NCNAXISp ) NC5=NC5+1
            END IF
         END IF
      END IF

                  END DO
               END DO
            END DO
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C CROSS PRODUCT
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CROSSPROD(Va,Vb,Vc)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Va(3),Vb(3),Vc(3)

      Vc(1) = Va(2)*Vb(3) -Va(3)*Vb(2)
      Vc(2) = Va(3)*Vb(1) -Va(1)*Vb(3)
      Vc(3) = Va(1)*Vb(2) -Va(2)*Vb(1)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C MAX COMMON DIVISOR (MCD)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FINDMCD(N,M,MCD)
      IMPLICIT REAL*8 (A-H,O-Z)

      MCD=1
      DO I=1,N
         IF ( MOD(N,I).EQ.0 .AND. MOD(M,I).EQ.0 ) MCD=I
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PRINT THE POINT GROUP FOUND
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE WRTPG(LUW,PG)
      IMPLICIT REAL*8 (A-H,O-Z)

      CHARACTER*4 PG

      WRITE (LUW,"(/A/A,A/A/)")
     !" -------------------------",
     !" * Point group found: ",PG,
     !" -------------------------"

      RETURN
      END  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C PRINT THE POINT GROUP FOUND
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NEWCN(V,RN,RK,NCNAXIS,CNAXIS)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DIMENSION V(3), CNAXIS(MAXOPSYM,5)

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0

      IF ( NCNAXIS.EQ.0 )  THEN
         NCNAXIS=1
         CNAXIS(NCNAXIS,1) = V(1)
         CNAXIS(NCNAXIS,2) = V(2)
         CNAXIS(NCNAXIS,3) = V(3)
c C_n^k : n = NMAX e k = ICN
         ICN=NINT(RK)
         NMAX=NINT(RN)
         CALL FINDMCD (NMAX,ICN,MCD)
         CNAXIS(NCNAXIS,4) = DFLOAT(NMAX/MCD)
         CNAXIS(NCNAXIS,5) = DFLOAT(ICN/MCD)
      ELSE
c checking if this new Cn is actually a "new" Cn
         NDIFF=0
         DO K=1,NCNAXIS
            DOT = V(1)*CNAXIS(K,1) + V(2)*CNAXIS(K,2) + V(3)*CNAXIS(K,3)
c More than 1 degree is considered not collinear
            IF (DABS(DOT).LT.DCOS(PI/180.D0)) NDIFF=NDIFF+1
c Axes can be collinear but their order are different
            THETA =2.D0*PI*RK/RN
            THETAp=2.D0*PI*CNAXIS(K,5)/CNAXIS(K,4)
            IF (DABS(DOT).GE.DCOS(PI/180.D0) .AND.
     ! DABS(THETAp-THETA).GT.(PI/180.D0) )  NDIFF=NDIFF+1
         END DO
         IF ( NDIFF.EQ.NCNAXIS ) THEN
            IF ( NCNAXIS.GT.MAXOPSYM ) THEN
               WRITE (*,'(//A,I4/A/A,I4//)')
     !" Error: The number of proper rotations found : ", NCNAXIS,
     ! "       is greater than the maximum order allowed for a ",
     ! "       point group ", MAXOPSYM
      STOP " PLACZEK will stops here "
            END IF
            NCNAXIS=NCNAXIS+1
            CNAXIS(NCNAXIS,1) = V(1)
            CNAXIS(NCNAXIS,2) = V(2)
            CNAXIS(NCNAXIS,3) = V(3)
            ICN=NINT(RK)
            NMAX=NINT(RN)
            CALL FINDMCD (NMAX,ICN,MCD)
            CNAXIS(NCNAXIS,4) = DFLOAT(NMAX/MCD)
            CNAXIS(NCNAXIS,5) = DFLOAT(ICN/MCD)
         END IF
      END IF

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MMATRIZ(A,B,C,NLA,MLA,NCA,MCA,NLB,MLB,NCB,MCB)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(MLA,MCA),B(MLB,MCB),C(MLA,MCB)

      IF (NCA.NE.NLB) THEN
          WRITE (*,'(/A,2I4,/A,2I4,/A/)')
     !'@ Matriz A possui dimensao ',NLA,NCA,
     !'@ Matriz B possui dimensao ',NLB,NCB,
     !'@ Portanto elas nao podem ser multiplicadas!'
      ENDIF

      DO I = 1,NLA
      DO J = 1,NCB
         C(I,J) = 0.D0
         DO K = 1,NCA
            C(I,J) = C(I,J) + A(I,K) * B(K,J)
         ENDDO
      ENDDO
      ENDDO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
C ROTATIONAL ANALYSIS
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC  
      SUBROUTINE ROTOANA (LUW)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION Tn(3,3), D(3), PMOM(3)

      BOHR    = 0.529177208D0
      PI      = 3.14159265358979323846D0
      PLANCK  = 6.626076D-34
      Avog    = 6.022137D+23
      CmLight = 2.99792458D+10

      LINEAR = .FALSE.

*     WRITE (LUW,'(/A/A/A/)')
*    !' *************************************',
*    !' **** OUTPUT FROM PROGRAM ROTOANA ****',
*    !' *************************************'
     
C STARTING EVALUATION OF INERTIA TENSOR
      DO I = 1,3
         DO J = 1,3
           Tn(I,J) = 0.D0
         END DO
      END DO

C EVALUATING CENTER OF MASS COORDINATES: CMX,CMY,CMZ
      CMX   = 0.D0
      CMY   = 0.D0
      CMZ   = 0.D0
      SMASS = 0.D0
      DO I = 1,NATOM
         IF (AU) THEN
            XYZ(I,1) = XYZ(I,1) * BOHR
            XYZ(I,2) = XYZ(I,2) * BOHR
            XYZ(I,3) = XYZ(I,3) * BOHR
         END IF
         AMASS = RMASS(I)
         CMX = CMX + AMASS * XYZ(I,1)
         CMY = CMY + AMASS * XYZ(I,2)
         CMZ = CMZ + AMASS * XYZ(I,3)
         SMASS = SMASS + AMASS
      END DO
      CMX = CMX / SMASS
      CMY = CMY / SMASS
      CMZ = CMZ / SMASS


C INERTIA TENSOR
      DO I = 1,NATOM
         AMASS = RMASS(I)
         CX = XYZ(I,1) - CMX
         CY = XYZ(I,2) - CMY
         CZ = XYZ(I,3) - CMZ

         Tn(1,1) = Tn(1,1) + AMASS * ( CY**2 + CZ**2 )
         Tn(1,2) = Tn(1,2) - AMASS *   CX*CY
         Tn(1,3) = Tn(1,3) - AMASS *   CX*CZ

         Tn(2,1) = Tn(1,2)
         Tn(2,2) = Tn(2,2) + AMASS * ( CX**2 + CZ**2 )
         Tn(2,3) = Tn(2,3) - AMASS *   CY*CZ

         Tn(3,1) = Tn(1,3)
         Tn(3,2) = Tn(2,3)
         Tn(3,3) = Tn(3,3) + AMASS * ( CX**2 + CY**2 )
      END DO


C DIAGONALIZATION OF INERTIA TENSOR

      WRITE (LUW,'(29X,A/A,3(F12.6,1X))')
     !'x            y            z',
     !' # Center of Mass (in A):',CMX,CMY,CMZ

      WRITE (LUW,'(/A/A)')
     !' Inertia Tensor (in u*A^2)',
     !' ------------------------- '
      WRITE (LUW,'(19X,"x",14X,"y",14X,"z")')
      WRITE (LUW,'(2X,"x",2X,3F15.6)')(Tn(1,J),J=1,3)
      WRITE (LUW,'(2X,"y",2X,3F15.6)')(Tn(2,J),J=1,3)
      WRITE (LUW,'(2X,"z",2X,3F15.6)')(Tn(3,J),J=1,3)

      CALL JACOB (Tn,D,3,3)

      DO I=1,3
         PMOM(I) = D(I)
      END DO
      CALL REORD (Tn,D,3,3)

C SETTING IA,IB,IC

      dIC = D(1)
      dIB = D(2)
      dIA = D(3)

C PRINTING IA,IB,IC

      WRITE (LUW,'(/A/A)')
     !' Principal Moments and Principal Axes of Inertia',
     !' -----------------------------------------------'
      WRITE (LUW,'(33X,"x       y       z")')
      WRITE (LUW,'(A,F15.6,A,3X,"(",3(F7.4,1X),")")')
     !' IA  =  ',dIA,'  u*A^2 ',(Tn(i,3),i=1,3)
      WRITE (LUW,'(A,F15.6,A,3X,"(",3(F7.4,1X),")")')
     !' IB  =  ',dIB,'  u*A^2 ',(Tn(i,2),i=1,3)
      WRITE (LUW,'(A,F15.6,A,3X,"(",3(F7.4,1X),")")')
     !' IC  =  ',dIC,'  u*A^2 ',(Tn(i,1),i=1,3)

C SEARCHING FOR LINEARITY
      Thres = 1.D-04
      WRITE (LUW,*)
      IF ( DABS(dIA).LT.Thres .AND. DABS(dIB-dIC).LT.Thres )
     !LINEAR = .TRUE.

C PRINTING ROTATIONAL CONSTANTS

      WRITE (LUW,'(/A/A)')
     !' Rotational Constants',
     !' --------------------'

      CONV = ( PLANCK * Avog ) / ( 1.D-23 * CmLight * 8.D0 * PI**2 )

      IF (.NOT.LINEAR) THEN

           B_IA = CONV / dIA
           B_IB = CONV / dIB
           B_IC = CONV / dIC
           B_IAhz = B_IA * CmLight / 1D+06
           B_IBhz = B_IB * CmLight / 1D+06
           B_IChz = B_IC * CmLight / 1D+06

           WRITE (LUW,10)B_IAhz,B_IA
           WRITE (LUW,11)B_IBhz,B_IB
           WRITE (LUW,12)B_IChz,B_IC

      ELSE

           B_IC = CONV / dIC
           B_IChz = B_IC * CmLight / 1D+06

           WRITE (LUW,11)B_IChz,B_IC
             
      END IF
      WRITE (LUW,*)

 10   FORMAT(1X,'A = ',1F15.2,' MHz     (',1F16.6,' cm-1)')
 11   FORMAT(1X,'B = ',1F15.2,' MHz     (',1F16.6,' cm-1)')
 12   FORMAT(1X,'C = ',1F15.2,' MHz     (',1F16.6,' cm-1)')

C ROTATIONAL SYMMETRY

      t1 = DABS( dIA - dIB )
      t2 = DABS( dIA - dIC )
      t3 = DABS( dIB - dIC )

      IF ( LINEAR .AND. NATOM.EQ.2 ) WRITE (LUW,'(A)')
     !' * Molecule Recognized as Linear (Diatomic)'

      IF ( LINEAR .AND. NATOM.GT.2 ) WRITE (LUW,'(A)')
     !' * Molecule Recognized as Linear '

      IF ( .NOT. LINEAR ) WRITE (LUW,'(A)')
     !' * Molecule Recognized as Non-Linear'

C ROTOR TYPE

      IF ( t1.LE.Thres .AND. t2.LE.Thres .AND. t3.LE.Thres
     !     .AND. .NOT.LINEAR )
     !WRITE (LUW,'(A)')' * The molecule is a Spherical Top'

      IF ( t1.GT.Thres .AND. t3.GT.Thres .AND. dIA.GT.Thres
     !     .AND. .NOT.LINEAR )
     !WRITE (LUW,'(A)')' * The molecule is an Asymmetric Top'

      IF ( ABS(dIC-dIB-dIA).LT.Thres .AND. .NOT.LINEAR )
     !WRITE (LUW,'(A)')' * The molecule is Planar'

      IF ( ABS(dIB-dIC).LT.Thres .AND. (dIB-dIA).GT.Thres
     !    .AND. dIC.GE.Thres .AND. .NOT.LINEAR ) 
     !WRITE (LUW,'(A)')' * The molecule is a Prolate Symmetric Top'

      IF ( ABS(dIA-dIB).LT.Thres .AND. (dIC-dIB).GT.Thres
     !    .AND. dIC.GE.Thres .AND. .NOT.LINEAR ) 
     !WRITE (LUW,'(A)')' * The molecule is a Oblate Symmetric Top'

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Matrix Multiplication 
C
      SUBROUTINE PMATRIX (A,B,C,N1,N2,N3) 
      IMPLICIT REAL*8 (A-H,O-Z)
              
      DIMENSION A(N1,N2), B(N2,N3), C(N1,N3) 

      DO ila = 1,N1
         DO icb = 1,N3
            C(ila, icb) = 0.D0
            DO i = 1,N2
               C(ila, icb) = C(ila,icb) + A(ila,i) * B(i,icb)
            END DO

         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C JACOB Matrix Diagonalization
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE JACOB(V,E,NDIM,INDEX)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION E(INDEX), D(INDEX,INDEX), V(INDEX,INDEX)

      CALL NEWJACOBI (V,D,NDIM,INDEX)

      DO I = 1,NDIM
         E(I) = V(I,I)
      END DO
      DO I = 1,NDIM
         DO J = 1,NDIM
            V(I,J) = D(I,J)
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C SUBROUTINE JACOB
C WRITEN BY Prof. Y. Hase
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE NEWJACOBI(H,D,N,INDEX)
      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4 (I-N)

      DIMENSION H(INDEX,INDEX),D(INDEX,INDEX)

      do i=1,n
         do j=1,n
            D(i,j)=0.d0
         end do
         D(i,i)=1.d0
      end do
      IF (N .EQ. 1) RETURN
      s=0.d0
      do i=1,n-1
         do j=i+1,n
            s=s+H(i,j)**2
         end do
      end do
      IF (s .lt. 1.D-20) RETURN
      t=dsqrt(s+s)
  110 t=t/DFLOAT(n)
  120 m=1
      do i=1,n-1
         do j=i+1,n
            if (dabs(H(i,j)).lt.t) goto 130
            m=2
            hii=H(i,i)
            hjj=H(j,j)
            hij=H(i,j)
            tang1=dsign(2.d0,(hii-hjj))*hij
            tang2=dabs(hii-hjj)+dsqrt((hii-hjj)**2+hij**2*4.d0)
            tang=tang1/tang2
            cosi=1.d0/dsqrt(1.d0+tang**2)
            sine=tang*cosi
            if (i.gt.1) then
               do k=1,i-1
                  hki=H(k,i)
                  H(k,i)= hki*cosi+H(k,j)*sine
                  H(k,j)=-hki*sine+H(k,j)*cosi
               end do
            end if
            H(i,i)=hii*cosi**2+hjj*sine**2+hij*sine*cosi*2.d0
            H(j,j)=hii*sine**2+hjj*cosi**2-hij*sine*cosi*2.d0
            if (i+1.le.j-1) then
               do k=i+1,j-1
                  hik=H(i,k)
                  H(i,k)= hik*cosi+H(k,j)*sine
                  H(k,j)=-hik*sine+H(k,j)*cosi
               end do
            end if
            H(i,j)=0.d0
            if (j.lt.n) then
               do k=j+1,n
                  hik=H(i,k)
                  H(i,k)= hik*cosi+H(j,k)*sine
                  H(j,k)=-hik*sine+H(j,k)*cosi
               end do
            end if
            do k=1,n
C               dki=D(i,k)
C               D(i,k)= dki*cosi+D(j,k)*sine
C               D(j,k)=-dki*sine+D(j,k)*cosi
               dki=D(k,i)
               D(k,i)= dki*cosi+D(k,j)*sine
               D(k,j)=-dki*sine+D(k,j)*cosi
            end do
  130    end do
      end do
      if (m.eq.2) goto 120
      IF (t .GT. 1.D-25) GO TO 110

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE REORD(V,D,NDIM,INDEX)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION V(INDEX,NDIM), D(INDEX), E(INDEX), VEC(INDEX,INDEX)

      DO I = 1,NDIM
         K = 0
         E(I) = -1.D+4
         DO J = 1,NDIM
            VEC(J,I) = 0.D0
            IF (E(I) .LT. D(J)) THEN
                E(I) = D(J)
                K = J
            END IF
         END DO
         D(K) = -1.D+5
         DO J = 1,NDIM
            VEC(J,I) = V(J,K)
         END DO
      END DO
      DO I = 1,NDIM
         D(I) = E(I)
         DO J = 1,NDIM
            V(I,J) = VEC(I,J)
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  CRC Handbook of Chemistry and Physics
C  Chief Editor:      David R. Lide  
C  Associate Editor:  H. P. R. Frederikse
C  78th edition; 1997-1998; CRC Press; USA
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION ABISTP(NAT,NISOTOP,IWAB)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION AM(1100)
      CHARACTER*(*) IWAB

      ABISTP = 0.D0

      DATA (AM(I),I=1,906) /
C        Z           A. MASS           ABUNDANCE      
     &   1.D0,       1.007825D0,     99.985000D0,   
     &   1.D0,       2.014102D0,      0.015000D0,   
     &   1.D0,       3.016049D0,     -1.000000D0,   
     &   2.D0,       4.002603D0,     99.999863D0,   
     &   2.D0,       3.016029D0,      0.000137D0,   
     &   3.D0,       7.016003D0,     92.500000D0,   
     &   3.D0,       6.015121D0,      7.500000D0,   
     &   4.D0,       9.012182D0,    100.000000D0,   
     &   5.D0,      11.009305D0,     80.100000D0,   
     &   5.D0,      10.012937D0,     19.900000D0,   
     &   6.D0,      12.000000D0,     98.900000D0,   
     &   6.D0,      13.003355D0,      1.100000D0,   
     &   6.D0,      14.003242D0,     -1.000000D0,   
     &   7.D0,      14.003074D0,     99.634000D0,   
     &   7.D0,      15.000109D0,      0.366000D0,   
     &   8.D0,      15.994915D0,     99.762000D0,   
     &   8.D0,      17.999160D0,      0.200000D0,   
     &   8.D0,      16.999131D0,      0.038000D0,   
     &   9.D0,      18.998403D0,    100.000000D0,   
     &  10.D0,      19.992436D0,     90.480000D0,   
     &  10.D0,      21.991383D0,      9.250000D0,   
     &  10.D0,      20.993843D0,      0.270000D0,   
     &  11.D0,      22.989770D0,    100.000000D0,   
     &  12.D0,      23.985042D0,     78.990000D0,   
     &  12.D0,      25.982594D0,     11.010000D0,   
     &  12.D0,      24.985837D0,     10.000000D0,   
     &  13.D0,      26.981539D0,    100.000000D0,   
     &  14.D0,      27.976927D0,     92.230000D0,   
     &  14.D0,      28.976495D0,      4.670000D0,   
     &  14.D0,      29.973771D0,      3.100000D0,   
     &  15.D0,      30.973762D0,    100.000000D0,   
     &  16.D0,      31.972071D0,     95.020000D0,   
     &  16.D0,      33.967867D0,      4.210000D0,   
     &  16.D0,      32.971459D0,      0.750000D0,   
     &  16.D0,      35.967081D0,      0.020000D0,   
     &  17.D0,      34.968853D0,     75.770000D0,   
     &  17.D0,      36.965903D0,     24.230000D0,   
     &  18.D0,      39.962384D0,     99.600000D0,   
     &  18.D0,      35.967546D0,      0.337000D0,   
     &  18.D0,      37.962733D0,      0.063000D0,   
     &  19.D0,      38.963707D0,     93.258100D0,   
     &  19.D0,      40.961825D0,      6.730200D0,   
     &  19.D0,      39.963999D0,      0.011700D0,   
     &  20.D0,      39.962591D0,     96.941000D0,   
     &  20.D0,      43.955481D0,      2.086000D0,   
     &  20.D0,      41.958618D0,      0.647000D0,   
     &  20.D0,      47.952533D0,      0.187000D0,   
     &  20.D0,      42.958767D0,      0.135000D0,   
     &  20.D0,      45.953693D0,      0.004000D0,   
     &  21.D0,      44.955910D0,    100.000000D0,   
     &  22.D0,      47.947947D0,     73.800000D0,   
     &  22.D0,      45.952629D0,      8.000000D0,   
     &  22.D0,      46.951764D0,      7.300000D0,   
     &  22.D0,      48.947871D0,      5.500000D0,   
     &  22.D0,      49.944792D0,      5.400000D0,   
     &  23.D0,      50.943962D0,     99.750000D0,   
     &  23.D0,      49.947161D0,      0.250000D0,   
     &  24.D0,      51.940510D0,     83.789000D0,   
     &  24.D0,      52.940651D0,      9.501000D0,   
     &  24.D0,      49.946046D0,      4.345000D0,   
     &  24.D0,      53.938883D0,      2.365000D0,   
     &  25.D0,      54.938047D0,    100.000000D0,   
     &  26.D0,      55.934939D0,     91.720000D0,   
     &  26.D0,      53.939613D0,      5.800000D0,   
     &  26.D0,      56.935396D0,      2.100000D0,   
     &  26.D0,      57.933277D0,      0.280000D0,   
     &  27.D0,      58.933198D0,    100.000000D0,   
     &  28.D0,      57.935346D0,     68.077000D0,   
     &  28.D0,      59.930788D0,     26.223000D0,   
     &  28.D0,      61.928346D0,      3.634000D0,   
     &  28.D0,      60.931058D0,      1.140000D0,   
     &  28.D0,      63.927968D0,      0.926000D0,   
     &  29.D0,      62.929599D0,     69.170000D0,   
     &  29.D0,      64.927793D0,     30.830000D0,   
     &  30.D0,      63.929145D0,     48.600000D0,   
     &  30.D0,      65.926035D0,     27.900000D0,   
     &  30.D0,      67.924846D0,     18.800000D0,   
     &  30.D0,      66.927129D0,      4.100000D0,   
     &  30.D0,      69.925325D0,      0.600000D0,   
     &  31.D0,      68.925580D0,     60.108000D0,   
     &  31.D0,      70.924701D0,     39.892000D0,   
     &  32.D0,      73.921177D0,     35.940000D0,   
     &  32.D0,      71.922080D0,     27.660000D0,   
     &  32.D0,      69.924250D0,     21.230000D0,   
     &  32.D0,      72.923463D0,      7.730000D0,   
     &  32.D0,      75.921402D0,      7.440000D0,   
     &  33.D0,      74.921594D0,    100.000000D0,   
     &  34.D0,      79.916520D0,     49.610000D0,   
     &  34.D0,      77.917308D0,     23.780000D0,   
     &  34.D0,      75.919212D0,      9.360000D0,   
     &  34.D0,      81.916698D0,      8.730000D0,   
     &  34.D0,      76.919913D0,      7.630000D0,   
     &  34.D0,      73.922475D0,      0.890000D0,   
     &  35.D0,      78.918336D0,     50.690000D0,   
     &  35.D0,      80.916289D0,     49.310000D0,   
     &  36.D0,      83.911507D0,     57.000000D0,   
     &  36.D0,      85.910616D0,     17.300000D0,   
     &  36.D0,      81.913482D0,     11.600000D0,   
     &  36.D0,      82.914135D0,     11.500000D0,   
     &  36.D0,      79.916380D0,      2.250000D0,   
     &  36.D0,      77.920396D0,      0.350000D0,   
     &  37.D0,      84.911794D0,     72.165000D0,   
     &  37.D0,      86.909187D0,     27.835000D0,   
     &  38.D0,      87.905619D0,     82.580000D0,   
     &  38.D0,      85.909267D0,      9.860000D0,   
     &  38.D0,      86.908884D0,      7.000000D0,   
     &  38.D0,      83.913430D0,      0.560000D0,   
     &  39.D0,      88.905849D0,    100.000000D0,   
     &  40.D0,      89.904703D0,     51.450000D0,   
     &  40.D0,      93.906315D0,     17.380000D0,   
     &  40.D0,      91.905039D0,     17.150000D0,   
     &  40.D0,      90.905644D0,     11.220000D0,   
     &  40.D0,      95.908275D0,      2.800000D0,   
     &  41.D0,      92.906377D0,    100.000000D0,   
     &  42.D0,      97.905407D0,     24.130000D0,   
     &  42.D0,      95.904679D0,     16.680000D0,   
     &  42.D0,      94.905841D0,     15.920000D0,   
     &  42.D0,      91.906808D0,     14.840000D0,   
     &  42.D0,      99.907477D0,      9.630000D0,   
     &  42.D0,      96.906021D0,      9.550000D0,   
     &  42.D0,      93.905085D0,      9.250000D0,   
     &  43.D0,      97.907215D0,     -1.000000D0,   
     &  44.D0,     101.904349D0,     31.600000D0,   
     &  44.D0,     103.905424D0,     18.700000D0,   
     &  44.D0,     100.905582D0,     17.000000D0,   
     &  44.D0,      98.905939D0,     12.700000D0,   
     &  44.D0,      99.904219D0,     12.600000D0,   
     &  44.D0,      95.907599D0,      5.520000D0,   
     &  44.D0,      97.905287D0,      1.880000D0,   
     &  45.D0,     102.905500D0,    100.000000D0,   
     &  46.D0,     105.903478D0,     27.330000D0,   
     &  46.D0,     107.903895D0,     26.460000D0,   
     &  46.D0,     104.905079D0,     22.330000D0,   
     &  46.D0,     109.905167D0,     11.720000D0,   
     &  46.D0,     103.904029D0,     11.140000D0,   
     &  46.D0,     101.905634D0,      1.020000D0,   
     &  47.D0,     106.905092D0,     51.839000D0,   
     &  47.D0,     108.904757D0,     48.161000D0,   
     &  48.D0,     113.903357D0,     28.730000D0,   
     &  48.D0,     111.902758D0,     24.130000D0,   
     &  48.D0,     110.904182D0,     12.800000D0,   
     &  48.D0,     109.903005D0,     12.490000D0,   
     &  48.D0,     112.904400D0,     12.220000D0,   
     &  48.D0,     115.904754D0,      7.490000D0,   
     &  48.D0,     105.906461D0,      1.250000D0,   
     &  48.D0,     107.904176D0,      0.890000D0,   
     &  49.D0,     114.903880D0,     95.700000D0,   
     &  49.D0,     112.904061D0,      4.300000D0,   
     &  50.D0,     119.902199D0,     32.590000D0,   
     &  50.D0,     117.901609D0,     24.230000D0,   
     &  50.D0,     115.901747D0,     14.530000D0,   
     &  50.D0,     118.903310D0,      8.590000D0,   
     &  50.D0,     116.902956D0,      7.680000D0,   
     &  50.D0,     123.905274D0,      5.790000D0,   
     &  50.D0,     121.903440D0,      4.630000D0,   
     &  50.D0,     111.904826D0,      0.970000D0,   
     &  50.D0,     113.902784D0,      0.650000D0,   
     &  50.D0,     114.903348D0,      0.340000D0,   
     &  51.D0,     120.903821D0,     57.360000D0,   
     &  51.D0,     122.904216D0,     42.640000D0,   
     &  52.D0,     129.906229D0,     33.800000D0,   
     &  52.D0,     127.904463D0,     31.690000D0,   
     &  52.D0,     125.903314D0,     18.950000D0,   
     &  52.D0,     124.904433D0,      7.139000D0,   
     &  52.D0,     123.902823D0,      4.816000D0,   
     &  52.D0,     121.903054D0,      2.603000D0,   
     &  52.D0,     122.904271D0,      0.908000D0,   
     &  52.D0,     119.904048D0,      0.096000D0,   
     &  53.D0,     126.904473D0,    100.000000D0,   
     &  54.D0,     131.904144D0,     26.900000D0,   
     &  54.D0,     128.904780D0,     26.400000D0,   
     &  54.D0,     130.905072D0,     21.200000D0,   
     &  54.D0,     133.905395D0,     10.400000D0,   
     &  54.D0,     135.907219D0,      8.900000D0,   
     &  54.D0,     129.903509D0,      4.100000D0,   
     &  54.D0,     127.903531D0,      1.910000D0,   
     &  54.D0,     123.905894D0,      0.100000D0,   
     &  54.D0,     125.904281D0,      0.090000D0,   
     &  55.D0,     132.905429D0,    100.000000D0,   
     &  56.D0,     137.905232D0,     71.700000D0,   
     &  56.D0,     136.905812D0,     11.230000D0,   
     &  56.D0,     135.904553D0,      7.854000D0,   
     &  56.D0,     134.905665D0,      6.592000D0,   
     &  56.D0,     133.904486D0,      2.417000D0,   
     &  56.D0,     131.905042D0,      0.101000D0,   
     &  56.D0,     129.906282D0,     -2.000000D0,   
     &  57.D0,     138.906347D0,     99.909800D0,   
     &  57.D0,     137.907105D0,      0.090200D0,   
     &  58.D0,     139.905433D0,     88.480000D0,   
     &  58.D0,     141.909241D0,     11.080000D0,   
     &  58.D0,     137.905985D0,      0.250000D0,   
     &  58.D0,     135.907140D0,      0.190000D0,   
     &  59.D0,     140.907647D0,    100.000000D0,   
     &  60.D0,     141.907719D0,     27.130000D0,   
     &  60.D0,     143.910083D0,     23.800000D0,   
     &  60.D0,     145.913113D0,     17.190000D0,   
     &  60.D0,     142.909810D0,     12.180000D0,   
     &  60.D0,     144.912570D0,      8.300000D0,   
     &  60.D0,     147.916889D0,      5.760000D0,   
     &  60.D0,     149.920887D0,      5.640000D0,   
     &  61.D0,     144.912743D0,     -1.000000D0,   
     &  62.D0,     151.919729D0,     26.700000D0,   
     &  62.D0,     153.922206D0,     22.700000D0,   
     &  62.D0,     146.914895D0,     15.000000D0,   
     &  62.D0,     148.917181D0,     13.800000D0,   
     &  62.D0,     147.914819D0,     11.300000D0,   
     &  62.D0,     149.917273D0,      7.400000D0,   
     &  62.D0,     143.911998D0,      3.100000D0,   
     &  63.D0,     152.921225D0,     52.200000D0,   
     &  63.D0,     150.919847D0,     47.800000D0,   
     &  64.D0,     157.924019D0,     24.840000D0,   
     &  64.D0,     159.927049D0,     21.860000D0,   
     &  64.D0,     155.922118D0,     20.470000D0,   
     &  64.D0,     156.923956D0,     15.650000D0,   
     &  64.D0,     154.922618D0,     14.800000D0,   
     &  64.D0,     153.920861D0,      2.180000D0,   
     &  64.D0,     151.919786D0,      0.200000D0,   
     &  65.D0,     158.925342D0,    100.000000D0,   
     &  66.D0,     163.929171D0,     28.200000D0,   
     &  66.D0,     161.926795D0,     25.500000D0,   
     &  66.D0,     162.928728D0,     24.900000D0,   
     &  66.D0,     160.926930D0,     18.900000D0,   
     &  66.D0,     159.925193D0,      2.340000D0,   
     &  66.D0,     157.924403D0,      0.100000D0,   
     &  66.D0,     155.924277D0,      0.060000D0,   
     &  67.D0,     164.930319D0,    100.000000D0,   
     &  68.D0,     165.930290D0,     33.600000D0,   
     &  68.D0,     167.932368D0,     26.800000D0,   
     &  68.D0,     166.932046D0,     22.950000D0,   
     &  68.D0,     169.935461D0,     14.900000D0,   
     &  68.D0,     163.929198D0,      1.610000D0,   
     &  68.D0,     161.928775D0,      0.140000D0,   
     &  69.D0,     168.934212D0,    100.000000D0,   
     &  70.D0,     173.938859D0,     31.800000D0,   
     &  70.D0,     171.936378D0,     21.900000D0,   
     &  70.D0,     172.938208D0,     16.120000D0,   
     &  70.D0,     170.936323D0,     14.300000D0,   
     &  70.D0,     175.942564D0,     12.700000D0,   
     &  70.D0,     169.934759D0,      3.050000D0,   
     &  70.D0,     167.933894D0,      0.130000D0,   
     &  71.D0,     174.940770D0,     97.410000D0,   
     &  71.D0,     175.942679D0,      2.590000D0,   
     &  72.D0,     179.946546D0,     35.100000D0,   
     &  72.D0,     177.943696D0,     27.297000D0,   
     &  72.D0,     176.943217D0,     18.606000D0,   
     &  72.D0,     178.945812D0,     13.629000D0,   
     &  72.D0,     175.941406D0,      5.206000D0,   
     &  72.D0,     173.940044D0,      0.162000D0,   
     &  73.D0,     180.947992D0,     99.988000D0,   
     &  73.D0,     179.947462D0,      0.012000D0,   
     &  74.D0,     183.950928D0,     30.670000D0,   
     &  74.D0,     185.954357D0,     28.600000D0,   
     &  74.D0,     181.948202D0,     26.300000D0,   
     &  74.D0,     182.950220D0,     14.300000D0,   
     &  74.D0,     179.946701D0,      0.130000D0,   
     &  75.D0,     186.955744D0,     62.600000D0,   
     &  75.D0,     184.952951D0,     37.400000D0,   
     &  76.D0,     191.961467D0,     41.000000D0,   
     &  76.D0,     189.958436D0,     26.400000D0,   
     &  76.D0,     188.958137D0,     16.100000D0,   
     &  76.D0,     187.955860D0,     13.300000D0,   
     &  76.D0,     186.955741D0,      1.600000D0,   
     &  76.D0,     185.953830D0,      1.580000D0,   
     &  76.D0,     183.952488D0,      0.020000D0,   
     &  77.D0,     192.962917D0,     62.700000D0,   
     &  77.D0,     190.960584D0,     37.300000D0,   
     &  78.D0,     194.964766D0,     33.800000D0,   
     &  78.D0,     193.962655D0,     32.900000D0,   
     &  78.D0,     195.964926D0,     25.300000D0,   
     &  78.D0,     197.967869D0,      7.200000D0,   
     &  78.D0,     191.961019D0,      0.790000D0,   
     &  78.D0,     189.959917D0,      0.010000D0,   
     &  79.D0,     196.966543D0,    100.000000D0,   
     &  80.D0,     201.970617D0,     29.860000D0,   
     &  80.D0,     199.968300D0,     23.100000D0,   
     &  80.D0,     198.968254D0,     16.870000D0,   
     &  80.D0,     200.970277D0,     13.180000D0,   
     &  80.D0,     197.966743D0,      9.970000D0,   
     &  80.D0,     203.973467D0,      6.870000D0,   
     &  80.D0,     195.965807D0,      0.150000D0,   
     &  81.D0,     204.974401D0,     70.476000D0,   
     &  81.D0,     202.972320D0,     29.524000D0,   
     &  82.D0,     207.976627D0,     52.400000D0,   
     &  82.D0,     205.974440D0,     24.100000D0,   
     &  82.D0,     206.975872D0,     22.100000D0,   
     &  82.D0,     203.973020D0,      1.400000D0,   
     &  83.D0,     208.980374D0,    100.000000D0,   
     &  84.D0,     208.982404D0,     -1.000000D0,   
     &  85.D0,     209.987126D0,     -1.000000D0,   
     &  86.D0,     222.017570D0,     -1.000000D0,   
     &  87.D0,     223.019733D0,     -1.000000D0,   
     &  88.D0,     226.025402D0,     -1.000000D0,   
     &  89.D0,     227.027750D0,     -1.000000D0,   
     &  90.D0,     232.038054D0,    100.000000D0,   
     &  91.D0,     231.035880D0,     -1.000000D0,   
     &  92.D0,     238.050785D0,     99.274500D0,   
     &  92.D0,     235.043924D0,      0.720000D0,   
     &  92.D0,     234.040947D0,      0.005500D0,   
     &  93.D0,     237.048168D0,     -1.000000D0,   
     &  94.D0,     239.052157D0,     -1.000000D0,   
     &  94.D0,     244.064199D0,     -1.000000D0,   
     & 119.D0,       1.D+06,        100.000000D0 /

      IST = 0
      DO I = 1,302
         J = 3 * I - 2
        IF ( NAT.EQ.NINT(AM(J)) ) IST = IST + 1
        IF ( IST.EQ.NISOTOP ) THEN
             IF ( IWAB.EQ.'MASS' )      ABISTP = AM(J+1)
             IF ( IWAB.EQ.'ABUNDANCE' ) ABISTP = AM(J+2)
             IF ( IWAB.NE.'MASS'.AND.IWAB.NE.'ABUNDANCE' )
     &            STOP 'ILLEGAL KEYWORD IN ABISTP'
             RETURN
        ENDIF
      ENDDO
      IF ( IST.GT.1 )  STOP 'NISOTOP exceeded in ABISTP'
      IF ( IST.EQ.0 ) STOP 'Z exceeded in ABISTP'

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Module to symmetrize molecules
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SYMMETRIZER(MAXCYCL)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION AVRDij(MAXAT,MAXAT)
      EXTERNAL func, dfunc
      COMMON /gradDelta/ AVRDij

      niter=0
      fretmax=1.d-4*DSQRT(DFLOAT(NATOM))
      itmax=2500
      fret=1.D3*fretmax
      DO WHILE ( fret.GT.fretmax .AND.niter.LT.itmax )
         CALL SMTRZR (iter,fret,MAXCYCL)
         niter=niter+iter
      END DO
      WRITE (80,"(/A,I6/A,1P,1D20.5)")
     *" Number of BFGS iterations used in symmetrization:",niter,
     *" Final value of 'delta(r_ij)' function: ",fret

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Computes THE INERTIA TENSOR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE ITENSOR(PIA,PIB,PIC,Tn)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION Tn(3,3), D(3)

c Calculating the Center of Mass 
      CMX = 0.D0
      CMY = 0.D0
      CMZ = 0.D0
      TMASS = 0.D0
      DO j = 1,NATOM
         AMASS = RMASS(j)
         CMX = CMX + XYZ(j,1) * AMASS
         CMY = CMY + XYZ(j,2) * AMASS
         CMZ = CMZ + XYZ(j,3) * AMASS
         TMASS = TMASS + AMASS
      END DO
      CMX = CMX / TMASS
      CMY = CMY / TMASS
      CMZ = CMZ / TMASS

C INERTIA TENSOR
      DO I = 1,3
         DO J = 1,3
           Tn(I,J) = 0.D0
         END DO
      END DO
      DO j = 1,NATOM
         CX = XYZ(j,1)-CMX
         CY = XYZ(j,2)-CMY
         CZ = XYZ(j,3)-CMZ 

         AMASS = RMASS(j)

         Tn(1,1) = Tn(1,1) + AMASS * ( CY**2 + CZ**2 )
         Tn(1,2) = Tn(1,2) - AMASS *   CX*CY
         Tn(1,3) = Tn(1,3) - AMASS *   CX*CZ
         Tn(2,1) = Tn(1,2)
         Tn(2,2) = Tn(2,2) + AMASS * ( CX**2 + CZ**2 )
         Tn(2,3) = Tn(2,3) - AMASS *   CY*CZ
         Tn(3,1) = Tn(1,3)
         Tn(3,2) = Tn(2,3)
         Tn(3,3) = Tn(3,3) + AMASS * ( CX**2 + CY**2 )
      END DO

C DIAGONALIZATION OF INERTIA TENSOR
      CALL JACOB (Tn,D,3,3)
      CALL REORD (Tn,D,3,3)

C SETTING IA,IB,IC
      PIC = D(1)
      PIB = D(2)
      PIA = D(3)

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Gradient of Delta({D_ij}) function
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE dfunc(p,GradDeltaRij)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON /XYZini/ XYZ0(MAXAT,3)

      DIMENSION DAiAj(MAXAT,MAXAT), AVRDij(MAXAT,MAXAT),
     *          GradDeltaRij(MAXNMOD), p(MAXNMOD)
      COMMON /gradDelta/ AVRDij

      DO i=1,3*NATOM+6
         GradDeltaRij(i) = 0.d0
      END DO

      k=0
      DO i=1,NATOM
         DO j=1,3
            k=k+1
            XYZ(i,j) = p(k)
         END DO
      END DO
    
c Interatomic distances (in Angstrom)
      DO I = 1,NATOM
         DO J = I,NATOM
            XX = ( XYZ(I,1) - XYZ(J,1) ) * ( XYZ(I,1) - XYZ(J,1) )
            YY = ( XYZ(I,2) - XYZ(J,2) ) * ( XYZ(I,2) - XYZ(J,2) )
            ZZ = ( XYZ(I,3) - XYZ(J,3) ) * ( XYZ(I,3) - XYZ(J,3) )
            DAiAj(I,J) = DSQRT( XX + YY + ZZ )
            DAiAj(J,I) = DAiAj(I,J)
         END DO
      END DO

      iGd=0
      DO i=1,NATOM
         DO iXYZ=1,3
            iGd=iGd+1
            DO j=1,NATOM
               IF ( j.NE.i ) THEN
       GradDeltaRij(iGd) = GradDeltaRij(iGd) + 
     * 2.d0*( DAiAj(i,j) - AVRDij(i,j) ) * 
     *               ( XYZ(i,iXYZ) - XYZ(j,iXYZ) ) / DAiAj(i,j) 
               END IF
       GradDeltaRij(iGd) = GradDeltaRij(iGd) + 
     * 2.d0*RMASS(j)*XYZ(j,iXYZ)*RMASS(i) 
               IF ( iXYZ.EQ.1 ) THEN
                  GradDeltaRij(iGd) = GradDeltaRij(iGd) + 
     * 2.d0*RMASS(j)*( XYZ0(j,1)*XYZ(j,2) -XYZ0(j,2)*XYZ(j,1) ) * 
     *      RMASS(i)*(          -XYZ0(i,2) )                    +
     * 2.d0*RMASS(j)*( XYZ0(j,1)*XYZ(j,3) -XYZ0(j,3)*XYZ(j,1) ) *
     *      RMASS(i)*(          -XYZ0(i,3) )
               ELSE IF ( iXYZ.EQ.2 ) THEN
                  GradDeltaRij(iGd) = GradDeltaRij(iGd) + 
     * 2.d0*RMASS(j)*( XYZ0(j,1)*XYZ(j,2) -XYZ0(j,2)*XYZ(j,1) ) * 
     *      RMASS(i)*(          +XYZ0(i,1) )                    +
     * 2.d0*RMASS(j)*( XYZ0(j,2)*XYZ(j,3) -XYZ0(j,3)*XYZ(j,2) ) *
     *      RMASS(i)*(          -XYZ0(i,3) )
               ELSE IF ( iXYZ.EQ.3 ) THEN
                  GradDeltaRij(iGd) = GradDeltaRij(iGd) + 
     * 2.d0*RMASS(j)*( XYZ0(j,1)*XYZ(j,3) -XYZ0(j,3)*XYZ(j,1) ) * 
     *      RMASS(i)*(          +XYZ0(i,1) )                    +
     * 2.d0*RMASS(j)*( XYZ0(j,2)*XYZ(j,3) -XYZ0(j,3)*XYZ(j,2) ) *
     *      RMASS(i)*(          +XYZ0(i,2) )
               END IF
            END DO
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C BFGS routine
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE dfpmin(p,n,gtol,iter,fret,func,dfunc)
      INTEGER iter,n,NNMAX,ITMAX, MAXAT

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, NNMAX=3*MAXAT)
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL*8 fret,gtol,p(NNMAX),func,EPS,STPMX,TOLX
      PARAMETER (ITMAX=500,STPMX=100.d0,EPS=3.D-8,TOLX=4.D0*EPS)
      EXTERNAL dfunc,func

C USES dfunc,func,lnsrch
c Given a starting point "p(1:n)" that is a vector of length "n", 
c the Broyden-Fletcher-Goldfarb-Shanno variant of Davidon-Fletcher-Powell 
c minimization is performed on a function "func", using its gradient as 
c calculated by a routine "dfunc". The convergence requirement on zeroing
c the gradient is input as "gtol". Returned quantities are : 
c "p(1:n)" (the location of the minimum),
c "iter" (the number of iterations that were performed), and
c "fret" (the minimum value of the function). The routine "lnsrch"
c is called to perform approximate line minimizations.
c
c Parameters:
c "NNMAX" is the maximum anticipated value of "n";
c "ITMAX" is the maximum allowed number of iterations;
c "STPMX" is the scaled maximum step length allowed in line searches;
c "TOLX" is the convergence criterion on "x" values.

      INTEGER i,its,j
      LOGICAL check
      REAL*8 den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp,test,
     *       dg(NNMAX),g(NNMAX),hdg(NNMAX),hessin(NNMAX,NNMAX),
     *       pnew(NNMAX),xi(NNMAX)
 
      fp=func(p)
c Calculate starting function value and gradient,
      call dfunc(p,g)
      sum=0.d0
      do i=1,n
c and initialize the inverse Hessian to the unit matrix.
         do j=1,n
            hessin(i,j)=0.d0
         enddo 
         hessin(i,i)=1.d0
         xi(i)=-g(i)
c Initial line direction.
         sum=sum+p(i)**2
      enddo 
      stpmax=STPMX*max(Dsqrt(sum),Dfloat(n))
      do its=1,ITMAX
c Main loop over the iterations.
         iter=its
         call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,func)
c The new function evaluation occurs in lnsrch; 
c save the function value in "fp" for the next
c line search. It is usually safe to ignore the value of check.
         fp=fret
         do i=1,n
            xi(i)=pnew(i)-p(i)
c Update the line direction,
            p(i)=pnew(i)
c and the current point.
         enddo 
         test=0.d0
c Test for convergence on "x".
         do i=1,n
            temp=Dabs(xi(i))/max(Dabs(p(i)),1.d0)
            if(temp.gt.test)test=temp
         enddo 
         if (test.lt.TOLX) return
         do i=1,n
c Save the old gradient,
            dg(i)=g(i)
         enddo 
         call dfunc(p,g)
c and get the new gradient.
         test=0.d0
c Test for convergence on zero gradient.
         den=max(fret,1.d0)
         do i=1,n
            temp=Dabs(g(i))*max(Dabs(p(i)),1.d0)/den
            if(temp.gt.test)test=temp
         enddo
         if (test.lt.gtol) return
         do i=1,n
c Compute diffence of gradients,
            dg(i)=g(i)-dg(i)
         enddo 
         do  i=1,n
c and diffence times current matrix.
            hdg(i)=0.d0
            do  j=1,n
               hdg(i)=hdg(i)+hessin(i,j)*dg(j)
            enddo 
         enddo 
         fac=0.d0
c Calculate dot products for the denominators.
         fae=0.d0
         sumdg=0.d0
         sumxi=0.d0
         do  i=1,n
            fac=fac+dg(i)*xi(i)
            fae=fae+dg(i)*hdg(i)
            sumdg=sumdg+dg(i)**2
            sumxi=sumxi+xi(i)**2
         enddo 
         if (fac.gt.Dsqrt(EPS*sumdg*sumxi)) then
c Skip update if "fac" not sumciently positive.
            fac=1.d0/fac
            fad=1.d0/fae
            do  i=1,n
c The vector that makes BFGS diffent from DFP:
               dg(i)=fac*xi(i)-fad*hdg(i)
            enddo
            do i=1,n
c The BFGS updating formula:
               do j=i,n
                  hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)
     * -fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j)
                  hessin(j,i)=hessin(i,j)
               enddo
            enddo 
         endif
         do i=1,n
c Now calculate the next direction to go,
            xi(i)=0.d0
            do j=1,n
               xi(i)=xi(i)-hessin(i,j)*g(j)
            enddo 
         enddo 
      enddo 
c and go back for another iteration.
c      pause 'too many iterations in dfpmin'
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE linmin(p,xi,n,fret)
      INTEGER n,NNMAX

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, NNMAX=3*MAXAT)
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL*8 fret,p(NNMAX),xi(NNMAX),TOL
      PARAMETER (TOL=1.d-4)
c Maximum anticipated "n" ,and "TOL" passed to brent.
C USES brent,f1dim,mnbrak 
c Given an n-dimensional point "p(1:n)" and an n-dimensional direction
c "xi(1:n)",moves and resets "p" to where the function "func(p)"
c takes on a minimum along the direction "xi" from "p", and replaces
c "xi" by the actual vector displacement that "p" was moved. Also returns c as "fret" the value of "func" at the returned location "p". 
c This is actually all accomplished by calling the routines "mnbrak"
c and "brent".
      INTEGER j,ncom
      REAL*8 ax,bx,fa,fb,fx,xmin,xx,pcom(NNMAX),xicom(NNMAX),brent,f1dim
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim
      ncom=n
c Set up the common block.
      do j=1,n
         pcom(j)=p(j)
         xicom(j)=xi(j)
      enddo
      ax=0.d0
c Initial guess for brackets.
      xx=1.d0
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do j=1,n
         xi(j)=xmin*xi(j)
         p(j)=p(j)+xi(j)
      enddo
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION f1dim(x)
      INTEGER NNMAX
      REAL*8 func,x

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, NNMAX=3*MAXAT)
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C USES func
c Used by "linmin" as the function passed to "mnbrak" and "brent".
      INTEGER j,ncom
      REAL*8 pcom(NNMAX),xicom(NNMAX),xt(NNMAX)
      COMMON /f1com/ pcom,xicom,ncom
      do j=1,ncom
         xt(j)=pcom(j)+x*xicom(j)
      enddo
      f1dim=func(xt)
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
      INTEGER n
      LOGICAL check

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, NNMAX=3*MAXAT)
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL*8 f,fold,stpmax,g(NNMAX),p(NNMAX),x(NNMAX),xold(NNMAX),func,
     *       ALF,TOLX
      PARAMETER (ALF=1.d-4,TOLX=1.d-7)
      EXTERNAL func
C USES func
c Given an n-dimensional point xold(1:n), the value of the function and 
c gradient there, fold and g(1:n) , and a direction p(1:n), 
c nds a new point x(1:n) along the direction p from xold where the 
c function func has decreased \sumciently." The new function value
c is returned in f.  
c stpmax is an input quantity that limits the length of the steps so that
c you do not try to evaluate the function in regions where it is 
c undened or subject to overflow.  p is usually the Newton direction. 
c The output quantity check is false on a normal exit.  It is true when x
c is too close to xold. In a minimization algorithm, this usually signals
c convergence and can be ignored. However, in a zero-nding algorithm thec calling program should check whether the convergence is spurious.
c Parameters: 
c ALF ensures sumcient decrease in function value;
c TOLX is the convergence criterion on x.
      INTEGER i
      REAL*8 a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,
     * sum,temp,test,tmplam
c vidal
      f2=0.d0
      alam2=0.d0
c vidal
      check=.false.
      sum=0.d0
      do i=1,n
         sum=sum+p(i)*p(i)
      enddo
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
c Scale if attempted step is too big.
        do i=1,n
           p(i)=p(i)*stpmax/sum
        enddo
      endif
      slope=0.d0
      do i=1,n
         slope=slope+g(i)*p(i)
      enddo
      if (slope.ge.0.d0) stop 'newsym.f: roundoff problem in lnsrch'
      test=0.d0
c Compute min.
      do i=1,n
         temp=Dabs(p(i))/max(Dabs(xold(i)),1.d0)
         if(temp.gt.test)test=temp
      enddo
      alamin=TOLX/test
      alam=1.d0
c Always try full Newton step rst.  1 continue
c Start of iteration loop.
 01   continue
      do i=1,n
         x(i)=xold(i)+alam*p(i)
      enddo
      f=func(x)
      if(alam.lt.alamin)then
c Convergence on x. For zero nding, the calling program should verify 
c the convergence.
         do i=1,n
            x(i)=xold(i)
         enddo
         check=.true.
         return
      else if(f.le.fold+ALF*alam*slope)then
c Sumficient function decrease.
         return
      else
c Backtrack.
         if(alam.eq.1.d0)then
c First time.
            tmplam=-slope/(2.d0*(f-fold-slope))
         else
c Subsequent backtracks.
            rhs1=f-fold-alam*slope
            rhs2=f2-fold-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/
     * (alam-alam2)
            if(a.eq.0.d0)then
              tmplam=-slope/(2.*b)
            else
              disc=b*b-3.d0*a*slope
              if(disc.lt.0.d0)then
                tmplam=.5d0*alam
              else if(b.le.0.d0)then
                 tmplam=(-b+dsqrt(disc))/(3.d0*a)
              else
                 tmplam=-slope/(b+dsqrt(disc))
              endif
            endif
            if(tmplam.gt..5d0*alam)tmplam=.5d0*alam
         endif
      endif
      alam2=alam
      f2=f
      alam=max(tmplam,.1d0*alam)
      goto 1
c Try again.
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      REAL*8 ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034d0, GLIMIT=100.d0, TINY=1.d-20)
c Given a function func, and given distinct initial points ax and bx, 
c this routine searches in the downhill direction (defined by the 
c function as evaluated at the initial points) and returns new points ax, 
c bx, cx that bracket a minimum of the function. Also returned are
cthe function values at the three points, fa, fb, and fc.
c Parameters:
c GOLD is the default ratio by which successive intervals are magnified;
c GLIMIT is the maximum magnication allowed for a parabolic-step.
      REAL*8 dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
cSwitch roles of a and b so that we can go downhill in the direction 
c from a to b.
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
cFirst guess for c.
      fc=func(cx)
 01   if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
c Compute u by parabolic extrapolation from a;b;c.  TINY is used to 
c prevent any possible division by zero.
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.d0*sign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
c We won't go farther than this. Test various possibilities:
        if((bx-u)*(u-cx).gt.0.d0)then
c Parabolic u is between b and c : try it.
           fu=func(u)
           if(fu.lt.fc)then
c Got a minimum between b and c.
             ax=bx
             fa=fb
             bx=u
             fb=fu
             return
           else if(fu.gt.fb)then
c Got a minimum between between a and u.
              cx=u
              fc=fu
              return
           endif
           u=cx+GOLD*(cx-bx)
c Parabolic t was no use. Use default magnication.
           fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.d0)then
c Parabolic t is between c and its allowed limit.
           fu=func(u)
           if(fu.lt.fc)then
             bx=cx
             cx=u
             u=cx+GOLD*(cx-bx)
             fb=fc
             fc=fu
             fu=func(u)
           endif
        else if((u-ulim)*(ulim-cx).ge.0.d0)then
c Limit parabolic u to maximum allowed value.
           u=ulim
           fu=func(u)
        else
c Reject parabolic u , use default magnication.
           u=cx+GOLD*(cx-bx)
           fu=func(u)
        endif
        ax=bx
c Eliminate oldest point and continue.
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      FUNCTION brent(ax,bx,cx,f,tol,xmin)
      INTEGER ITMAX
      REAL*8 brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=100,CGOLD=.3819660d0,ZEPS=1.0d-10)
c Given a function f, and given a bracketing triplet of abscissas ax,
c bx, cx (such that bx is between ax and cx ,and f(bx) is less than both 
c f(ax) and f(cx)), this routine isolates the minimum to a fractional 
c precision of about tol using Brent's method. The abscissa of the 
c minimum is returned as xmin , and the minimum function value is 
c returned as brent , the returned function value.  
c Parameters: 
c Maximum allowed number of iterations; golden ratio; and a small number 
c that protects against trying to achieve fractional accuracy for a 
c minimum that happens to be exactly zero.
      INTEGER iter
      REAL*8 a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
c vidal
      d=0.d0
c vidal
      a=min(ax,cx) 
c a and b must be in ascending order, though the input abscissas need 
c not be.
      b=max(ax,cx)
      v=bx
c Initializations...
      w=v
      x=v
      e=0.d0
c This will be the distance moved on the step before last.
      fx=f(x)
      fv=fx
      fw=fx
      do iter=1,ITMAX
c Main program loop.
        xm=0.5d0*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
c Test for done here.
        if(dabs(e).gt.tol1) then
c Construct a trial parabolic t.
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if(q.gt.0.d0) p=-p
          q=dabs(q)
          etemp=e
          e=d
          if(abs(p).ge.dabs(.5d0*q*etemp).or.p.le.q*(a-x).or.
     * p.ge.q*(b-x)) goto 1
c The above conditions determine the acceptability of the parabolic t. 
c Here it is o.k.:
            d=p/q
c Take the parabolic step.
            u=x+d
            if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
            goto 2
c Skip over the golden section step.
          endif
 1        if(x.ge.xm) then
c We arrive here for a golden section step, which we take
c into the larger of the two segments.
            e=a-x
          else
            e=b-x
          endif
          d=CGOLD*e
c Take the golden section step.
 2        if(abs(d).ge.tol1) then
c Arrive here with d computed either from parabolic t, or
c else from golden section.
            u=x+d
          else
            u=x+sign(tol1,d)
          endif
          fu=f(u)
c This is the one function evaluation per iteration,
          if(fu.le.fx) then
c and now we have to decide what to do with our function
c evaluation. Housekeeping follows:
            if(u.ge.x) then
              a=x
            else
              b=x
            endif
            v=w
            fv=fw
            w=x
            fw=fx
            x=u
            fx=fu
          else
            if(u.lt.x) then
              a=u
            else
              b=u
            endif
            if(fu.le.fw .or. w.eq.x) then
              v=w
              fv=fw
              w=u
              fw=fu
            else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
              v=u
              fv=fu
            endif
          endif
c Done with housekeeping. Back for another iteration.
      enddo
      stop 'newsym.f: brent exceed maximum iterations'
 3    xmin=x
c Arrive here ready to exit with best values.  
      brent=fx
      return
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Quadratic deviation function: 
C Delta({D_ij}) = \sum_{i=1,j>i}^{NATOM} ( D_ij - <D_ij> )**2 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DOUBLE PRECISION FUNCTION func(p)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON /XYZini/ XYZ0(MAXAT,3)

      DIMENSION Dij(MAXAT,MAXAT), AVRDij(MAXAT,MAXAT),p(MAXNMOD)
      COMMON /gradDelta/ AVRDij

      n=0
      do i=1,natom
         do j=1,3
            n=n+1
            XYZ(i,j) = p(n)
         end do
      end do

c Interatomic distances 
      DO I = 1,NATOM
         DO J = I,NATOM
            XX = ( XYZ(I,1) - XYZ(J,1) ) * ( XYZ(I,1) - XYZ(J,1) )
            YY = ( XYZ(I,2) - XYZ(J,2) ) * ( XYZ(I,2) - XYZ(J,2) )
            ZZ = ( XYZ(I,3) - XYZ(J,3) ) * ( XYZ(I,3) - XYZ(J,3) )
            Dij(I,J) = DSQRT( XX + YY + ZZ )
            Dij(J,I) = Dij(I,J)
         END DO
      END DO

      DeltaRij =0.d0
      DO i = 1,NATOM-1
         DO j = i+1,NATOM
            DeltaRij = DeltaRij + (Dij(i,j)-AVRDij(i,j))**2
         END DO
      END DO
      CMX=0.D0
      CMY=0.D0
      CMZ=0.D0
      dIxy=0.d0
      dIxz=0.d0
      dIyz=0.d0
      DO j = 1,NATOM
         CMX = CMX + RMASS(j) * XYZ(j,1)
         CMY = CMY + RMASS(j) * XYZ(j,2)
         CMZ = CMZ + RMASS(j) * XYZ(j,3)
         dIxy = dIxy + 
     *   RMASS(j) * (XYZ0(j,1)*XYZ(j,2) -XYZ0(j,2)*XYZ(j,1))
         dIxz = dIxz + 
     *   RMASS(j) * (XYZ0(j,1)*XYZ(j,3) -XYZ0(j,3)*XYZ(j,1))
         dIyz = dIyz + 
     *   RMASS(j) * (XYZ0(j,2)*XYZ(j,3) -XYZ0(j,3)*XYZ(j,2))
      END DO
      func=DeltaRij +  CMX**2 +  CMY**2 +  CMZ**2 +
     *                dIxy**2 + dIxz**2 + dIyz**2

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c Grouping atoms into sets of Symmetrically Equivalent Ones (SEA)
c SU : symmetry-unique atom
c SE : symmetrically equivalent atom
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE FINDSEA(DAiAj,Thres,EQLDij,LINELWST)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION DAiAj(MAXAT,MAXAT), EQLDij(MAXAT,MAXAT), 
     *          Di(MAXAT), Dj(MAXAT),LINELWST(MAXAT,MAXAT)
      LOGICAL   SEA(MAXAT), JUMP(MAXAT)

      DO i=1,MAXAT
         DO j=1,MAXAT
            ISEA(i,j) = 0
            LINELWST(i,j) = 0
            EQLDij(i,j) = 0
         END DO
         SEA(i)  = .FALSE.
         NSEA(i) = 0
      END DO
      nseaset = 0
      RATOM = DFLOAT(NATOM)

      DO iSU = 1,NATOM

         IF ( .NOT.SEA(iSU) ) THEN
            nseaset = nseaset + 1
            NSEA(nseaset) = 1
            ISEA(nseaset,NSEA(nseaset)) = iSU
            SEA(iSU) = .TRUE.
            DO linSU=1,NATOM
               EQLDij(linSU,iSU) = DAiAj(linSU,iSU) 
               Di(linSU)         = DAiAj(linSU,iSU)
            END DO
         END IF
         RMASSiSU=RMASS(iSU)

         iSE=iSU+1
         DO WHILE ( iSE.LE.NATOM ) 
            NEQD = 0
            DO linSE=1,NATOM
               Dj(linSE) = DAiAj(linSE,iSE)
            END DO
            RNORMDIFF = DABS( DiNORM(Di,NATOM,MAXAT) -
     *                        DiNORM(Dj,NATOM,MAXAT) ) / RATOM
            IF ( RNORMDIFF.LE.Thres .AND. RMASSiSU.EQ.RMASS(iSE) 
     *           .AND. .NOT.SEA(iSE) ) THEN

               DO linSE=1,NATOM
                  JUMP(linSE) = .FALSE.
               END DO
               DO linSU=1,NATOM
                  RLWST = 1.d10
                  LWSTlinSE =  0
                  RMASSlinSU=RMASS(linSU)
                  DilinSU=Di(linSU)
                  DO linSE=1,NATOM
                     DIFFline=DABS(DilinSU-Dj(linSE))
                     IF ( .NOT. JUMP(linSE) 
     *                    .AND. DIFFline  .LE.Thres 
     *                    .AND. RMASSlinSU.EQ.RMASS(linSE)
     *                    .AND. DIFFline  .LT.RLWST         
     *                  ) THEN
                        RLWST = DIFFline
                        LWSTlinSE = linSE
                     END IF
                  END DO
                  IF ( RLWST.LT.1.d10 ) THEN
                     NEQD = NEQD+1 
                     JUMP(LWSTlinSE) = .TRUE.
                     EQLDij(linSU,iSE) = Dj(LWSTlinSE) 
                     LINELWST(linSU,iSE) = LWSTlinSE
                     Dj(LWSTlinSE) = -1.D0
                  END IF
               END DO

            END IF

            IF ( NEQD.EQ.NATOM ) THEN
               NSEA(nseaset) = NSEA(nseaset) + 1
               ISEA(nseaset,NSEA(nseaset)) = iSE
               SEA(iSE) = .TRUE.
            END IF
            iSE=iSE+1
         END DO

      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C This routine computes the average interatomic distance matrix <D_ij>
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE AVERAGEdij(Thres,EQLDij,LINELWST,AVRDij,IERROR)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      DIMENSION EQLDij(MAXAT,MAXAT), 
     *          LINELWST(MAXAT,MAXAT), 
     *          AVRDij(MAXAT,MAXAT), tmpDij(MAXAT,MAXAT)
      LOGICAL   JUMP(MAXAT)

c Ordering sets of SEA according to their number of atoms 
c going from the more populated to the less.
c      DO i=1,nseaset
c         JUMP(i) = .FALSE.
c      END DO
c      DO i=1,nseaset
c         ISORTSEA(i) = 0
c         LGST = 0
c         DO J=1,nseaset
c            IF ( NSEA(J).GT.LGST .AND. .NOT.JUMP(J) ) THEN
c               LGST = NSEA(J)
c               ISORTSEA(i) = J
c            END IF
c         END DO
c         JUMP(ISORTSEA(i)) = .TRUE.
c      END DO

c Computing the averaged interatomic distances AVRD_ij = <D_ij> 
c
c LINELWST : line where the lowest value for the difference 
c            D_iSU - D_iSEA were found
c D_iSU : interatomic distance associated to a symmetry-unique atom
c D_iSE : interatomic distance of a symmetrically-equivalent atom
      ZEROrij=0.7D0
      IF (AU) ZEROrij=ZEROrij/0.529177208D0
      DO i=1,nseaset
         iSU = ISEA(i,1)
         DO linSU=1,NATOM
            AVRDij(linSU,iSU) = 0.d0
            IDEL=0
            DO k=1,NSEA(i)
               iSE = ISEA(i,k)
               IF      ( EQLDij(linSU,iSE).GT.ZEROrij ) THEN
                  AVRDij(linSU,iSU) = AVRDij(linSU,iSU) + 
     *                                EQLDij(linSU,iSE)
               ELSE IF ( EQLDij(linSU,iSE).LT.ZEROrij .AND. 
     *                   linSU.NE.iSU                 )  THEN
                  IDEL=IDEL+1
               END IF
            END DO
            IF (IDEL.LT.NSEA(i) ) THEN
               AVRDij(linSU,iSU) = AVRDij(linSU,iSU) / 
     *                             DFLOAT(NSEA(i)-IDEL)
            END IF
            IF ( AVRDij(linSU,iSU).LT.ZEROrij .AND. linSU.NE.iSU ) THEN
               WRITE (*,"(//A,I4,A,I4,A,F10.6/A,F7.3,A/)") 
     *" Warning: Average distance between atoms",linSU," and ",iSU,
     *" is too small:",AVRDij(linSU,iSU),
     *" Note: Values below ",ZEROrij," angstrom are considered zero."
               IERROR=1
               RETURN
            END IF
            JUMP(linSU) = .FALSE.
         END DO

c Averaging over distances appearing in lines of Dij that are  
c interpreted as equal (those lower than "Thres" and that corresponds 
c to a pair of atoms having same mass) but in fact are slightly 
c different.
c
c e.g.:
c       O  H1  H2
c O     0
c H1   d1
c H2   d2
c
c d1 = (d1+d2) / 2 ; d2 = d1
c
         DO linSU=1,NATOM
            JUMP(linSU)=.FALSE.
            DO icolSU=1,NATOM
               tmpDij(linSU,icolSU) = AVRDij(linSU,iSU)
            END DO
         END DO
         DO linSU=1,NATOM-1
            NeqDij = 0
            SumDij = 0.d0
            DO linSU2=linSU+1,NATOM
               IF ( DABS(tmpDij(linSU,iSU)-tmpDij(linSU2,iSU)).LT.Thres 
     !              .AND.  RMASS(linSU).EQ.RMASS(linSU2) 
     !              .AND. .NOT.JUMP(linSU2)  ) THEN 
                  SumDij = SumDij + AVRDij(linSU2,iSU)
                  NeqDij=NeqDij+1
                  JUMP(linSU2) = .TRUE.
               END IF
            END DO
            IF (NeqDij.GT.0) AVRDij(linSU,iSU) = SumDij / DFLOAT(NeqDij)
         END DO
         DO linSU=1,NATOM
            JUMP(linSU)=.FALSE.
         END DO
         DO linSU=1,NATOM-1
            DO linSU2=linSU+1,NATOM
               IF ( DABS(tmpDij(linSU,iSU)-tmpDij(linSU2,iSU)).LT.Thres 
     !              .AND.  RMASS(linSU).EQ.RMASS(linSU2) 
     !              .AND.  .NOT.JUMP(linSU2)  ) THEN 
                  JUMP(linSU2) = .TRUE.
                  AVRDij(linSU2,iSU) = AVRDij(linSU,iSU)
               END IF
            END DO
         END DO

c Filling columns of D_ij corresponding to SEA with the
c average values <r_ij> computed for the SU atom.
         IF ( NSEA(i).GT.1 ) THEN
            iSU = ISEA(i,1)
            DO k=2,NSEA(i)
               iSE = ISEA(i,k)
               DO linSU=1,NATOM
                  AVRDij(LINELWST(linSU,iSE),iSE) = AVRDij(linSU,iSU)
               END DO
            END DO
         END IF
      END DO

c The average distance (<r_ij> or <r_ji>) calculated using the 
c interatomic distances from the most populated SEA set 
c (that having more atoms) will replace the other one.
c      DO i=1,nseaset
c         incr = ISORTSEA(nseaset+1-i)
c         imax = ISORTSEA(1)
c         DO j=1,NSEA(incr)
c            iSE = ISEA(incr,j)
c            DO linSE=1,NATOM
c               IF ( incr.LT.imax) THEN
c
c                  LineAtm=0
c                  NatmSEA=0
c                  do ii=1,nseaset
c                  do jj=1,NSEA(ii)
c                     LineAtm=LineAtm+1
c                     IF ( LineAtm.EQ.linSE ) NatmSEA=NSEA(ii)
c                  end do
c                  end do
c                  IF ( NatmSEA.GT.NSEA(incr) ) 
c     *               AVRDij(linSE,iSE) = AVRDij(iSE,linSE)
c               END IF
c            END DO
c         END DO
c      END DO
c
c Now Dij is symmetrized
c      DO i=1,NATOM
c      DO j=i,NATOM
c         avr = ( AVRDij(i,j) + AVRDij(j,i) ) / 2.D0
c         AVRDij(i,j) = avr
c         AVRDij(j,i) = avr
c      END DO
c      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Routine grouping the main actions of symmetrizer 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SMTRZR(iter,fret,MAXCYCL)
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON /XYZini/ XYZ0(MAXAT,3)

      DIMENSION DAiAj(MAXAT,MAXAT), EQLDij(MAXAT,MAXAT),
     *          LINELWST(MAXAT,MAXAT), 
     *          AVRDij(MAXAT,MAXAT), Tn(3,3),
     *          XYZR(MAXAT,3),p(NNMAX),g(NNMAX),CM(3)
      EXTERNAL func, dfunc
      COMMON /gradDelta/ AVRDij

      BOHR = 0.529177208D0

c Calculating the Center of Mass 
      DO i=1,3
         CM(i) = 0.D0
      END DO
      TMASS = 0.D0
      DO I = 1,NATOM
         AMASS = RMASS(I)
         DO J=1,3
            CM(J) = CM(J) + XYZ(I,J) * AMASS
         END DO
         TMASS = TMASS + AMASS
      END DO
      DO i=1,3
         CM(i) = CM(i) / TMASS
      END DO

c Putting C.M. in the origin
      CONV = 1.D0
      IF (AU) CONV = BOHR
      DO I = 1,NATOM
         DO J=1,3
             XYZ(I,J) = ( XYZ(I,J)-CM(J) ) * CONV
         END DO
      END DO

c Aligning Principal Axes with X, Y and Z ones
      CALL ITENSOR (PIA,PIB,PIC,Tn)
      DO I3=1,NATOM
         DO J3=1,3
            XYZR(I3,J3) = 0.D0
            DO K3=1,3
               XYZR(I3,J3) = XYZR(I3,J3) + XYZ(I3,K3) * Tn(K3,J3)
            END DO
         END DO
      END DO
      DO I3=1,NATOM
         DO J3=1,3
            XYZ(I3,J3) = XYZR(I3,J3)
c            IF ( niter.EQ.0 ) XYZ0(I3,J3) = XYZR(I3,J3)
         END DO
      END DO

 01   CONTINUE
c Interatomic distances 
      DO I = 1,NATOM
         DO J = I,NATOM
            XX = ( XYZ(I,1) - XYZ(J,1) ) * ( XYZ(I,1) - XYZ(J,1) )
            YY = ( XYZ(I,2) - XYZ(J,2) ) * ( XYZ(I,2) - XYZ(J,2) )
            ZZ = ( XYZ(I,3) - XYZ(J,3) ) * ( XYZ(I,3) - XYZ(J,3) )
            DAiAj(I,J) = DSQRT( XX + YY + ZZ )
            DAiAj(J,I) = DAiAj(I,J)
         END DO
      END DO

c Grouping atoms into sets of Symmetrically Equivalent Ones (SEA)
c Here, several thresholds are tried in order to find the lowest
c number of SEA sets.
c We mean that a molecule having only one SEA set is more symmetric
c than another that has two or more sets. 
      MINSETS=NATOM+1
      ThresNew=1.D-3
      NCYCL=0
c      DO WHILE ( MINSETS.GT.NATOM .OR. NCYCL.LT.MAXCYCL )
      DO WHILE ( NCYCL.LT.MAXCYCL )
         NCYCL=NCYCL+1
         ThresNew = ThresNew * 1.8d0
         CALL FINDSEA (DAiAj,ThresNew,EQLDij,LINELWST)
         IF ( MINSETS.GT.nseaset ) THEN
            MINSETS = nseaset 
            Thres = ThresNew
         END IF
      END DO
      CALL FINDSEA (DAiAj,Thres,EQLDij,LINELWST)
      IF (AU) THEN
         WRITE (*,"(/A,1X,1F10.4,1X,A/A/)") 
     *" > Distances lower than",Thres,"bohr are considered equal.",
     *" Note: Only 'symmetrizer' module uses this threshold."
      ELSE
         WRITE (*,"(/A,1X,1F10.4,1X,A/A/)") 
     *" > Distances lower than",Thres,"angstrom are considered equal.",
     *" Note: Only 'symmetrizer' module uses this threshold."
      END IF

c Evaluating the averaged interatomic distance matrix: AVRDij = <r_ij>
      IERROR=0
      CALL AVERAGEdij (Thres,EQLDij,LINELWST,AVRDij,IERROR)
      IF ( IERROR.EQ.1 .OR. nseaset.EQ.NATOM ) THEN
         DO I3=1,NATOM
            DO J3=1,3
               XYZ(I3,J3) = XYZ0(I3,J3)
            END DO
         END DO
         MAXCYCL=MAXCYCL-1
         IF ( MAXCYCL.GT.5) GO TO 01
      END IF

c BFGS minimization
      n=0
      DO i=1,NATOM
         DO j=1,3
            n=n+1
            p(n) = XYZ(i,j)
         END DO
      END DO
      gtol = 1.d-10
      fini=func(p)
      CALL dfunc (p,g)
      CALL dfpmin (p,n,gtol,iter,fret,func,dfunc)
      n=0
      DO i=1,NATOM
         DO j=1,3
            n=n+1
            XYZ(i,j) = p(n)
         END DO
      END DO

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Starting coordinates
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE XYZzero
      IMPLICIT REAL*8 (A-H,O-Z)

c BLOCK SHARED WITH MANY ROUTINES cccccccccccccccccccccccccccccccccccccc
      PARAMETER (MAXAT=1000, MAXNMOD=3*MAXAT, NNMAX=MAXNMOD,
     *           MAXOPSYM=240)
      CHARACTER ATMSBL(MAXAT)*2
      LOGICAL   AU, LINEAR
      COMMON /SYMMETRY/ XYZ(MAXAT,3), RMASS(MAXAT), CHARGE(MAXAT),
     *                  ISEA(MAXAT,MAXAT), NSEA(MAXAT), NSEASET,
     *                  NATOM, ATMSBL, AU, LINEAR
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      COMMON /XYZini/ XYZ0(MAXAT,3)

      DIMENSION Tn(3,3), XYZR(MAXAT,3), CM(3)

      BOHR = 0.529177208D0

c Calculating the Center of Mass 
      DO i=1,3
         CM(i) = 0.D0
      END DO
      TMASS = 0.D0
      DO I = 1,NATOM
         AMASS = RMASS(I)
         DO J=1,3
            CM(J) = CM(J) + XYZ(I,J) * AMASS
         END DO
         TMASS = TMASS + AMASS
      END DO
      DO i=1,3
         CM(i) = CM(i) / TMASS
      END DO

c Putting C.M. in the origin
      CONV = 1.D0
      IF (AU) CONV = BOHR
      DO I = 1,NATOM
         DO J=1,3
             XYZ(I,J) = ( XYZ(I,J)-CM(J) ) * CONV
         END DO
      END DO

c Aligning Principal Axes with X, Y and Z ones
      CALL ITENSOR (PIA,PIB,PIC,Tn)
      DO I3=1,NATOM
         DO J3=1,3
            XYZR(I3,J3) = 0.D0
            DO K3=1,3
               XYZR(I3,J3) = XYZR(I3,J3) + XYZ(I3,K3) * Tn(K3,J3)
            END DO
         END DO
      END DO
      DO I3=1,NATOM
         DO J3=1,3
            XYZ0(I3,J3) = XYZR(I3,J3)
         END DO
      END DO

      RETURN
      END

