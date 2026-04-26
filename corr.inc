c
c......corr.inc - Parameters read or computed from the WFN file.
c
c  ndets  = Maximum number of Slater dets in multideterminantal WFNs.
c  ncore  = Number of core orbitals.
c  nact   = Number of valence or active orbitals.
c  nelact = Number of valence or active electrons.
c  nel    = Number of electrons.
c  nalpha = Number of alpha electrons.
c  nbeta  = Number of beta  electrons.
c  mocore = Number of core orbitals excluded from the calculation.
c  moval  = Number of valence orbitals included in the calculation.
c      default mocore=0, default moval=nmo (total number of MOs)
c
       common /icor/ ndets,ncore,nact,nelact,nel,nalpha,nbeta,
     &               mocore,moval
