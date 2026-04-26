c
c-----------------------------------------------------------------------
c
      SUBROUTINE nlmfill ()
      include  'param.inc'
c
c.....nlm keeps the nlm values of x^n y^l z^m gaussian
c
      do j=1,3
         do i=1,maxtype
            nlm(i,j)=0
         enddo
      enddo
c
c-----warning------warning------warning------warning------warning-------
c
c     WFN files coming from gamess02 must use the definitions of the 
c     NLM's currently commented. However, WFN files coming from gamess10 
c     or later versions must use the NON COMMENTED definitions.
c
c-----warning------warning------warning------warning------warning-------
c
c.....p's
c
C     nlm(2,1)=1      !px       gamess02
C     nlm(3,2)=1      !py       gamess02
C     nlm(4,3)=1      !pz       gamess02
c                     !         gamess02
c.....d's             !         gamess02
c                     !         gamess02
C     nlm(5,1)=2      !xx       gamess02
C     nlm(6,2)=2      !yy       gamess02
C     nlm(7,3)=2      !zz       gamess02
C     nlm(8,1)=1      !xy       gamess02
C     nlm(8,2)=1      !         gamess02
C     nlm(9,1)=1      !xz       gamess02
C     nlm(9,3)=1      !         gamess02
C     nlm(10,2)=1     !yz       gamess02
C     nlm(10,3)=1     !         gamess02
c                     !         gamess02
c.....f's             !         gamess02
c                     !         gamess02
C     nlm(11,1)=3     !xxx      gamess02
C     nlm(12,2)=3     !yyy      gamess02
C     nlm(13,3)=3     !zzz      gamess02
C     nlm(14,1)=2     !xxy      gamess02
C     nlm(14,2)=1     !         gamess02
C     nlm(15,1)=2     !xxz      gamess02
C     nlm(15,3)=1     !         gamess02
C     nlm(16,2)=2     !yyz      gamess02
C     nlm(16,3)=1     !         gamess02
C     nlm(17,1)=1     !xyy      gamess02
C     nlm(17,2)=2     !         gamess02
C     nlm(18,1)=1     !xzz      gamess02
C     nlm(18,3)=2     !         gamess02
C     nlm(19,2)=1     !yzz      gamess02
C     nlm(19,3)=2     !         gamess02
C     nlm(20,1)=1     !xyz      gamess02
C     nlm(20,2)=1     !         gamess02
C     nlm(20,3)=1     !         gamess02
c                     !         gamess02
c.....g's             !         gamess02
c                     !         gamess02
C     nlm(21,3)=4     !zzzz     gamess02
C     nlm(22,2)=1     !yzzz     gamess02
C     nlm(22,3)=3     !         gamess02
C     nlm(23,2)=2     !yyzz     gamess02
C     nlm(23,3)=2     !         gamess02
C     nlm(24,2)=3     !yyyz     gamess02
C     nlm(24,3)=1     !         gamess02
C     nlm(25,2)=4     !yyyy     gamess02
C     nlm(26,1)=1     !xzzz     gamess02
C     nlm(26,3)=3     !         gamess02
C     nlm(27,1)=1     !xyzz     gamess02
C     nlm(27,2)=1     !         gamess02
C     nlm(27,3)=2     !         gamess02
C     nlm(28,1)=1     !xyyz     gamess02
C     nlm(28,2)=2     !         gamess02
C     nlm(28,3)=1     !         gamess02
C     nlm(29,1)=1     !xyyy     gamess02
C     nlm(29,2)=3     !         gamess02
C     nlm(30,1)=2     !xxzz     gamess02
C     nlm(30,3)=2     !         gamess02
C     nlm(31,1)=2     !xxyz     gamess02
C     nlm(31,2)=1     !         gamess02
C     nlm(31,3)=1     !         gamess02
C     nlm(32,1)=2     !xxyy     gamess02
C     nlm(32,2)=2     !         gamess02
C     nlm(33,1)=3     !xxxz     gamess02
C     nlm(33,3)=1     !         gamess02
C     nlm(34,1)=3     !xxxy     gamess02
C     nlm(34,2)=1     !         gamess02
C     nlm(35,1)=4     !xxxx     gamess02
c
c
c.....p's
c
      nlm(2,1)=1      !px       gamess10
      nlm(3,2)=1      !py       gamess10
      nlm(4,3)=1      !pz       gamess10
c                     !         gamess10
c.....d's             !         gamess10
c                     !         gamess10
      nlm(5,1)=2      !xx       gamess10
      nlm(6,2)=2      !yy       gamess10
      nlm(7,3)=2      !zz       gamess10
      nlm(8,1)=1      !xy       gamess10
      nlm(8,2)=1      !         gamess10
      nlm(9,1)=1      !xz       gamess10
      nlm(9,3)=1      !         gamess10
      nlm(10,2)=1     !yz       gamess10
      nlm(10,3)=1     !         gamess10
c                     !         gamess10
c.....f's             !         gamess10
c                     !         gamess10
      nlm(11,1)=3     !xxx      gamess10
      nlm(12,2)=3     !yyy      gamess10
      nlm(13,3)=3     !zzz      gamess10
      nlm(14,1)=2     !xxy      gamess10
      nlm(14,2)=1     !         gamess10
      nlm(15,1)=2     !xxz      gamess10
      nlm(15,3)=1     !         gamess10
      nlm(16,2)=2     !yyz      gamess10
      nlm(16,3)=1     !         gamess10
      nlm(17,1)=1     !xyy      gamess10
      nlm(17,2)=2     !         gamess10
      nlm(18,1)=1     !xzz      gamess10
      nlm(18,3)=2     !         gamess10
      nlm(19,2)=1     !yzz      gamess10
      nlm(19,3)=2     !         gamess10
      nlm(20,1)=1     !xyz      gamess10
      nlm(20,2)=1     !         gamess10
      nlm(20,3)=1     !         gamess10
c                     !         gamess10
c.....g's             !         gamess10
c                     !         gamess10
      nlm(23,3)=4     !zzzz     gamess10
      nlm(29,2)=1     !yzzz     gamess10
      nlm(29,3)=3     !         gamess10
      nlm(32,2)=2     !yyzz     gamess10
      nlm(32,3)=2     !         gamess10
      nlm(27,2)=3     !yyyz     gamess10
      nlm(27,3)=1     !         gamess10
      nlm(22,2)=4     !yyyy     gamess10
      nlm(28,1)=1     !xzzz     gamess10
      nlm(28,3)=3     !         gamess10
      nlm(35,1)=1     !xyzz     gamess10
      nlm(35,2)=1     !         gamess10
      nlm(35,3)=2     !         gamess10
      nlm(34,1)=1     !xyyz     gamess10
      nlm(34,2)=2     !         gamess10
      nlm(34,3)=1     !         gamess10
      nlm(26,1)=1     !xyyy     gamess10
      nlm(26,2)=3     !         gamess10
      nlm(31,1)=2     !xxzz     gamess10
      nlm(31,3)=2     !         gamess10
      nlm(33,1)=2     !xxyz     gamess10
      nlm(33,2)=1     !         gamess10
      nlm(33,3)=1     !         gamess10
      nlm(30,1)=2     !xxyy     gamess10
      nlm(30,2)=2     !         gamess10
      nlm(25,1)=3     !xxxz     gamess10
      nlm(25,3)=1     !         gamess10
      nlm(24,1)=3     !xxxy     gamess10
      nlm(24,2)=1     !         gamess10
      nlm(21,1)=4     !xxxx     gamess10
      return
      end
