c
c-----------------------------------------------------------------------
c
c---------Compute and write entropy (S): 
c                                               |    p(na,nb,...)   |
c         S = SUM (na,nb,...) p(na,nb,...)  log | ----------------- |
c                                               |   p(na) p(nb)...  |
c--------------------------------------------------------------------
c
      subroutine mutent (pnew,resnc,mocogrp,mocogrpa,mocogrpb,
     &   nprobt,np,ngroup,nel,wrt,lw)
      include   'implicit.inc'
      real(kind=8), dimension (:)     :: pnew(nprobt)
      integer(kind=4), dimension (:,:):: resnc(nprobt,ngroup)
      integer(kind=4)  mocogrp(ngroup)
      integer(kind=4)  mocogrpa(ngroup)
      integer(kind=4)  mocogrpb(ngroup)
      parameter (mxg = 10)
      integer(kind=4), allocatable,dimension (:)       :: maxelec
      integer(kind=4), allocatable,dimension (:)       :: minelec
      real(kind=8),    allocatable,dimension (:,:)     :: pgroup
      integer(kind=4), allocatable,dimension (:,:)     :: indp
      real(kind=8),    allocatable,dimension (:,:)     :: pgt
      integer(kind=4), allocatable,dimension (:)       :: jj,igroup
      logical wrt,allgtz,alleq
      real(kind=8),    allocatable,dimension (:)       :: s1
c
c-----Begin computation of Mutual Entropy Information.
c
      allocate (maxelec(ngroup))
      allocate (minelec(ngroup))
      allocate (pgroup(ngroup,0:nel))
      allocate (indp(np,mxg))
      allocate (jj(mxg))
      allocate (igroup(mxg))
      allocate (s1(mxg))

      write (lw,677)
      minelec=+1000
      maxelec=-1000
      pgroup=0d0
c
c-----Determine probabilities of each group p_1(na)
c
      do n=1,np
        do igr=1,ngroup
          ncurr=resnc(n,igr)+mocogrpa(igr)+mocogrpb(igr)
          pgroup(igr,ncurr)=pnew(n)+pgroup(igr,ncurr)
          minelec(igr)=min(minelec(igr),ncurr)
          maxelec(igr)=max(maxelec(igr),ncurr)
        enddo
      enddo
c
c-----Write probabilities of each group p_1(n_a)
c
      write (lw,678) 'MinElec = ',(minelec(igr),igr=1,ngroup)
      write (lw,678) 'MaxElec = ',(maxelec(igr),igr=1,ngroup)
      s1=0d0
      do igr=1,ngroup
        write (lw,*)'# Probabilities > 1D-8 for group ',igr
        do i=minelec(igr),maxelec(igr)
          probal=pgroup(igr,i)
          adds1=0d0
          if (probal.gt.0d0) adds1=probal*log(probal)
          s1(igr)=s1(igr)-adds1
          if (probal.ge.1D-8) write (lw,683) i,probal
        enddo
        write (lw,*)'# Shannon entropy  = ',s1(igr)
      enddo
c
c-----There are more than a group
c
      if (ngroup.gt.1) then
        ngprod=ngroup*(ngroup-1)/2
        allocate (pgt(ngprod,np))
        pgt=0d0
        i=0
        do igr1=1,ngroup
          igroup(1)=igr1
          do igr2=igr1+1,ngroup
            igroup(2)=igr2
            i=i+1
            write (lw,684) igr1,igr2
            ns=0
            do n=1,np
              if (pnew(n).ge.0d0) then
                jj(1:2)=resnc(n,igroup(1:2)) + 
     &            mocogrpa(igroup(1:2))+mocogrpb(igroup(1:2))
                if (ns.eq.0) then
                  ns = 1
                  pgt(i,ns) = pnew(n)
                  indp(ns,1:2) = jj(1:2)
                else
                  do m=1,ns
                    alleq=.true.
                    do l=1,2
                      alleq=alleq.and.jj(l).eq.indp(m,l)
                    enddo
                    if (alleq) then
                      pgt(i,m) = pgt(i,m)+pnew(n)
                      goto 1
                    endif
                  enddo
                  ns = ns+1
                  pgt(i,ns) = pnew(n)
                  indp(ns,1:2) = jj(1:2)
 1                continue
                endif
              endif
            enddo
            entropy=0d0
            st=0d0
            do m=1,ns
              jj(1:2) = indp(m,1:2)
              ptotal  = pgt(i,m)
              pall = 1d0
              allgtz = .true.
              do l=1,2
                pall = pall*pgroup(igroup(l),jj(l))
                allgtz = allgtz.and.pgroup(igroup(l),jj(l)).gt.0d0
              enddo
              addpt=0d0
              if (ptotal.gt.0d0) addpt=ptotal*log(ptotal)
              st=st-addpt
              if (ptotal.gt.0d0.and.allgtz) then
                entropy = entropy + ptotal*log(ptotal/pall)
                if (ptotal.gt.1D-6) then
                  write (lw,681) (jj(l),l=1,2),ptotal,pall,
     &              (pgroup(igroup(l),jj(l)),l=1,2)
                endif
              endif
            enddo
            write (lw,33) entropy,st
          enddo
        enddo
        deallocate (pgt)
      endif
c
c-----There are more than two group
c
      if (ngroup.gt.2) then
        ngprod=ngroup*(ngroup-1)*(ngroup-2)/6
        allocate (pgt(ngprod,np))
        pgt=0d0
        i=0
        do igr1=1,ngroup
          igroup(1)=igr1
          do igr2=igr1+1,ngroup
            igroup(2)=igr2
            do igr3=igr2+1,ngroup
              igroup(3)=igr3
              i=i+1
              write (lw,684) igr1,igr2,igr3
              ns=0
              do n=1,np
                if (pnew(n).ge.0d0) then
                  jj(1:3)=resnc(n,igroup(1:3)) +
     &               mocogrpa(igroup(1:3))+mocogrpb(igroup(1:3))
                  if (ns.eq.0) then
                    ns = 1
                    pgt(i,ns)   = pnew(n)
                    indp(ns,1:3) = jj(1:3)
                  else
                    do m=1,ns
                      alleq=.true.
                      do l=1,3
                        alleq=alleq.and.jj(l).eq.indp(m,l)
                      enddo
                      if (alleq) then
                        pgt(i,m) = pgt(i,m)+pnew(n)
                        goto 2
                      endif
                    enddo
                    ns = ns+1
                    pgt(i,ns)   = pnew(n)
                    indp(ns,1:3) = jj(1:3)
 2                  continue
                  endif
                endif
              enddo
              entropy=0d0
              st=0d0
              do m=1,ns
                jj(1:3) = indp(m,1:3)
                ptotal  = pgt(i,m)
                pall = 1d0
                allgtz = .true.
                do l=1,3
                  pall = pall*pgroup(igroup(l),jj(l))
                  allgtz = allgtz.and.pgroup(igroup(l),jj(l)).gt.0d0
                enddo
                addpt=0d0
                if (ptotal.gt.0d0) addpt=ptotal*log(ptotal)
                st=st-addpt
                if (ptotal.gt.0d0.and.allgtz) then
                  entropy = entropy + ptotal*log(ptotal/pall)
                  if (ptotal.gt.1D-6) then
                    write (lw,6810) (jj(l),l=1,3),ptotal,pall,
     &                (pgroup(igroup(l),jj(l)),l=1,3)
                  endif
                endif
              enddo
              write (lw,33) entropy,st
            enddo
          enddo
        enddo
        deallocate (pgt)
      endif
c
c-----There are more than three group
c
      if (ngroup.gt.3) then
        ngprod=ngroup*(ngroup-1)*(ngroup-2)*(ngroup-3)/24
        allocate (pgt(ngprod,np))
        pgt=0d0
        i=0
        do igr1=1,ngroup
          igroup(1)=igr1
          do igr2=igr1+1,ngroup
            igroup(2)=igr2
            do igr3=igr2+1,ngroup
              igroup(3)=igr3
              do igr4=igr3+1,ngroup
                igroup(4)=igr4
                i=i+1
                write (lw,684) igr1,igr2,igr3,igr4
                ns=0
                do n=1,np
                  if (pnew(n).ge.0d0) then
                    jj(1:4)=resnc(n,igroup(1:4))+
     &                 mocogrpa(igroup(1:4))+mocogrpb(igroup(1:4))
                    if (ns.eq.0) then
                      ns = 1
                      pgt(i,ns)   = pnew(n)
                      indp(ns,1:4) = jj(1:4)
                    else
                      do m=1,ns
                        alleq=.true.
                        do l=1,4
                          alleq=alleq.and.jj(l).eq.indp(m,l)
                        enddo
                        if (alleq) then
                          pgt(i,m) = pgt(i,m)+pnew(n)
                          goto 3
                        endif
                      enddo
                      ns = ns+1
                      pgt(i,ns)   = pnew(n)
                      indp(ns,1:4) = jj(1:4)
 3                  continue
                    endif
                  endif
                enddo
                entropy=0d0
                st=0d0
                do m=1,ns
                  jj(1:4) = indp(m,1:4)
                  ptotal  = pgt(i,m)
                  pall = 1d0
                  allgtz = .true.
                  do l=1,4
                    pall = pall*pgroup(igroup(l),jj(l))
                    allgtz = allgtz.and.pgroup(igroup(l),jj(l)).gt.0d0
                  enddo
                  addpt=0d0
                  if (ptotal.gt.0d0) addpt=ptotal*log(ptotal)
                  st=st-addpt
                  if (ptotal.gt.0d0.and.allgtz) then
                    entropy = entropy + ptotal*log(ptotal/pall)
                    if (ptotal.gt.1D-6) then
                      write (lw,6811) (jj(l),l=1,4),ptotal,pall,
     &                  (pgroup(igroup(l),jj(l)),l=1,4)
                    endif
                  endif
                enddo
                write (lw,33) entropy,st
              enddo
            enddo
          enddo
        enddo
        deallocate (pgt)
      endif
c
c-----There are more than three group
c
      if (ngroup.gt.4) then
        ngprod=ngroup*(ngroup-1)*(ngroup-2)*(ngroup-3)*(ngroup-4)/120
        allocate (pgt(ngprod,np))
        pgt=0d0
        i=0
        do igr1=1,ngroup
          igroup(1)=igr1
          do igr2=igr1+1,ngroup
            igroup(2)=igr2
            do igr3=igr2+1,ngroup
              igroup(3)=igr3
              do igr4=igr3+1,ngroup
                igroup(4)=igr4
                do igr5=igr4+1,ngroup
                  igroup(5)=igr5
                  i=i+1
                  write (lw,684) igr1,igr2,igr3,igr4,igr5
                  ns=0
                  do n=1,np
                    if (pnew(n).ge.0d0) then
                      jj(1:5)=resnc(n,igroup(1:5))+
     &                    mocogrpa(igroup(1:5))+mocogrpb(igroup(1:5))
                      if (ns.eq.0) then
                        ns = 1
                        pgt(i,ns) = pnew(n)
                        indp(ns,1:5) = jj(1:5)
                      else
                        do m=1,ns
                          alleq=.true.
                          do l=1,5
                            alleq=alleq.and.jj(l).eq.indp(m,l)
                          enddo
                          if (alleq) then
                            pgt(i,m) = pgt(i,m)+pnew(n)
                            goto 4
                          endif
                        enddo
                        ns = ns+1
                        pgt(i,ns)   = pnew(n)
                        indp(ns,1:5) = jj(1:5)
 4                    continue
                      endif
                    endif
                  enddo
                  entropy=0d0
                  st=0d0
                  do m=1,ns
                    jj(1:5) = indp(m,1:5)
                    ptotal  = pgt(i,m)
                    pall = 1d0
                    allgtz = .true.
                    do l=1,5
                      pall = pall*pgroup(igroup(l),jj(l))
                      allgtz = allgtz.and.pgroup(igroup(l),jj(l)).gt.0d0
                    enddo
                    addpt=0d0
                    if (ptotal.gt.0d0) addpt=ptotal*log(ptotal)
                    st=st-addpt
                    if (ptotal.gt.0d0.and.allgtz) then
                      entropy = entropy + ptotal*log(ptotal/pall)
                      if (ptotal.gt.1D-6) then
                        write (lw,6812) (jj(l),l=1,5),ptotal,pall,
     &                    (pgroup(igroup(l),jj(l)),l=1,5)
                      endif
                    endif
                  enddo
                  write (lw,33) entropy,st
                enddo
              enddo
            enddo
          enddo
        enddo
        deallocate (pgt)
      endif
      deallocate (maxelec,minelec,pgroup,indp,jj,igroup)
c
c-----End computation of Mutual Entropy Information.
c                                                   
      return
 677  format (/,' # ',64('-'),
     &  /,' # Mutual information entropy between groups',/,' # ',
     & 64('-'))
 678  format (' # ',a,2x,20I4)
 681  format (' # p(',2I3,' ) = ',F9.6,' ?= ',F9.6,' = ',F9.6,
     & ' * ', F9.6)
 6810 format (' # p(',3I3,' ) = ',F9.6,' ?= ',F9.6,' = ',F9.6,
     & ' * ', F9.6,' * ',F9.6)
 6811 format (' # p(',4I3,' ) = ',F9.6,' ?= ',F9.6,' = ',F9.6,
     & ' * ', F9.6,' * ',F9.6,' * ', F9.6)
 6812 format (' # p(',5I3,' ) = ',F9.6,' ?= ',F9.6,' = ',F9.6,
     & ' * ', F9.6,' * ',F9.6,' * ', F9.6,' * ', F9.6)
 33   format (' # Mutual Information Entropy = ',E15.8,/,
     &        ' # Shannon entropy            = ',E15.8)
 683  format (1x,'# p(',I3,') = ',F25.18)
 684  format (/1x,'# Probabilities > 1D-6 for groups ',10(I4,' +'))
      end
