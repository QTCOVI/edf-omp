c
c-----------------------------------------------------------------------
c   _ ________   ____
c  | |___|___  \| ___|
c  | |__ | | \  | |__     
c  |  __|| |  | |  __|
c  | |___|_|_/  | |
c  |_|___|_____/|_|
c
c-----------------------------------------------------------------------
c
      program edf
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! |  ______    ______    ______     ______    ______    ______     ______    ______    ______         |
! | /_____/\  /_____/\  /_____/\   /_____/\  /_____/\  /_____/\   /_____/\  /_____/\  /_____/\        |
! | \::::_\/_ \:::_ \ \ \::::_\/_  \::::_\/_ \:::_ \ \ \::::_\/_  \::::_\/_ \:::_ \ \ \::::_\/_       |
! |  \:\/___/\ \:\ \ \ \ \:\/___/\  \:\/___/\ \:\ \ \ \ \:\/___/\  \:\/___/\ \:\ \ \ \ \:\/___/\      |
! |   \::___\/  \:\ \ \ \ \:::._\/   \::___\/  \:\ \ \ \ \:::._\/   \::___\/  \:\ \ \ \ \:::._\/      |
! |    \:\____/\ \:\/.:| | \:\ \      \:\____/\ \:\/.:| | \:\ \      \:\____/\ \:\/.:| | \:\ \        |
! |     \_____\/  \____/_/  \_\/       \_____\/  \____/_/  \_\/       \_____\/  \____/_/  \_\/        |
! |                                                                                                   |
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ELECTRON NUMBER DISTRIBUTION FUNCTIONS
! (c) E. Francisco & A. Martín Pendas, 2022
!       University of Oviedo
! 
c     THIRD VERSION OF THE 'edf' PROGRAM. 
c
c     The following tasks are carried out:
c
c     1) Read the wavefunction and the atomic overlap matrix (AOM) 
c     2) Compute the electron number distribution function (EDF)
c        using the fragments defined in the input file.
c     3) Localize the input natural orbitals using the Cioslowski 
c        isopycnic localization method.
c     4) Compute the EDF from the wave function expanded in terms
c        of the isopycnic MOs. If all works properly, it should be 
c        exactly the same than that obtained in 2)
c     5) Compute an approximate EDF assuming that all the localized 
c        orbitals are orthogonal within all the defined fragments. This 
c        is exact when there are only two fragments and the localized 
c        orbitals are those resulting from the diagonalization of the 
c        overlap matrix (within one of the two fragments) between the 
c        molecular orbitals corresponding to a single-determinant wave 
c        function.
c
c        Tasks 3,4, and 5 are optional.
c
c.....EDF reads the aimpac .WFN file written by the GAMESS or GAUSSIAN 
c     packages, reads (or computes) the overlap integrals between the 
c     molecular orbitals (MOs) within all the basins, and obtains the 
c     probabilities of having exactly an integer number of electrons
c     n_1 in fragment 1, n_2 electrons in  fragment 2, ..., and n_m 
c     electrons in fragment 'm' of the molecule. The sum of the 'm' 
c     fragments must be necessarily equal to the whole molecule, and 
c     each atom can only be ascribed to a single fragment.
c
c.....Three categories of basins can be used in EDF:
c
c     1) The atomic basins of the Quantum Theory of Atoms in Molecules
c        (QTAIM) or atomic QTAIM basins.
c     2) The basins of the Electron Localization Function (ELF) or
c        ELF basins.
c     3) The basins (computed within the program itself) corresponding
c        to a partition of the molecule into fuzzy atomic basins. There
c        are three of such fuzzy atomic partitions included in the EDF 
c        code which are based on Mulliken, Lowdin, and MinDef atoms, 
c        respectively. Mulliken and Lowdin atoms are well-known and do
c        not deserve further comments. 'MinDef' atoms correspond to the
c        Minimally Deformed Atoms in Molecules defined by J. Fernandez 
c        rico et. al.. See for instance, 
c
c           a) "Analysis of the molecular density" 
c               Fernández Rico, J.; López, R.; Ramírez, G. 
c               J Chem Phys  1999, 110, 4213-4220.
c
c           b) "Analysis of the molecular density: STO densities2
c               Fernández Rico, J.; López, R.; Ema, I.; Ramírez, G.  
c               J Chem Phys 2002, 117, 533-540. 
c
c.....The overlap integrals within the atomic QTAIM basins are read in
c     by EDF with three different formats, corresponding to three 
c     different computational codes: PROMOLDEN, PROAIM, and TOPMOD.
c
c.....The overlap integrals within the ELF basins are read in exclusi-
c     vely with the format of TOPMOD code.
c
c
c.....Input (standard unit).............................................
c     (CAPITAL lettters are fixed names in the input)
c
c=====Record 1: ioverlap
c
c     ioverlap =   0  ==> The file with the overlaps between MOs comes from
c                         the PROMOLDEN program, or 'filedat' in the next 
c                         record is 'mulliken', 'lowdin','mindef', 'mindefrho',
c                         'netrho', 'promrho', 'heselmann', or 'becke'.
c
c     ioverlap =   1  ==> The file with the overlaps between MOs comes from 
c                         the PROAIM program.
c
c     ioverlap =   2  ==> The file with the overlaps between MOs comes from
c                         the TOPMOD program with ELF basins.
c
c     ioverlap =   3  ==> The file with the overlaps between MOs comes from
c                         the TOPMOD program with QTAIM basins.
c
c     ioverlap =   4  ==> The file with the overlaps between MOs comes from
c                         the AIMALL program with QTAIM basins.
c
c     ioverlap =  -1 ==>  Complex Atomic Overlap Matrix (See 'edfx.f' file
c                         for information about the input in this case)
c
c     ioverlap =  -2 ==>  Complex Atomic Overlap Matrix (See 'edfcrit.f' file
c                         for information about the input in this case)
c
c     ioverlap = -11 ==>  This option can be used when the wavefunction is
c                         a single determinant. After reading the AOM 
c                         integrals in the second record of the input file
c                         the control is transferred to the 'edfstd.f' 
c                         routine, where additional input can be read in
c                         (see the 'edfstd.f' file for this).
c
c=====Record 2: filedat  
c
c        Name of the file containing the atomic overlap matrix (AOM) 
c        between canonical MOs.
c
c        Instead of an existing 'filedat' file, the following keywords
c        are also valid in record 2.
c
c        filedat = 'mulliken' ---> The AOM is computed using a space 
c                   partitioning according to the Mulliken population 
c                   analysis method.
c        
c
c        filedat = 'lowdin' ---> The AOM is computed using a space 
c                   partitioning according to the Lowdin population 
c                   analysis method.
c
c        filedat = 'mindef' ---> The AOM is computed using a space 
c                   partitioning according to the minimally deformed 
c                   (mindef) atoms defined by Fernadez Rico et at.
c
c        filedat = 'mindefrho' ---> Same as above but with weights
c                   'w_A=rho_A/rho', that always give S_ii^A > 0 values.
c
c                   In this case, three integers can be given after 
c                   the MINDEFRHO keyword: NANG, NRAD and IRMESH,
c                   that represent the number of angular (NANG) and 
c                   radial (NRAD)  points in the integration, and the
c                   radial mapping (IRMESH).
c
c        filedat = 'netrho' ---> Same as above but with weights
c                   'w_A=rho_A/rho' determined with the so-called net
c                   densities, which elliminate all mixed (i.e. primitives
c                   in two different centers) densities when the w_A are
c                   computed.
c
c                   In this case, three integers can be given after 
c                   the MINDEFRHO keyword: NANG, NRAD and IRMESH,
c                   that represent the number of angular (NANG) and 
c                   radial (NRAD)  points in the integration, and the
c                   radial mapping (IRMESH).
c
c        filedat = 'promrho' ---> Same as above with 'w_A=rho_A/rho',
c                   being each 'rho_A' the IN VACUO atomic density 
c                   of the atom computed with the seven term analytical
c                   expression derived by fitting the Koga STO 
c                   atomic densities with the AUXPI code.
c
c        filedat = 'heselmann' ---> Identifies that the weights used
c                  in a chemical localization are those defined by
c                  Andreas Heselmann, JCTC 12, 2720-2741 (2016).
c
c        filedat = 'becke' ---> The AOM is computed using Becke's recipe
c
c                   In this case, four integers can be given
c                   after the MINDEFRHO keyword: NANG, NRAD, IRMESH, and 
c                   POW, that represent the number of angular (NANG) and 
c                   radial (NRAD)  points in the integration, the radial 
c                   mapping (IRMESH), and the iterative K parameter of
c                   Becke's method (POW).
c
c
c        In these five cases, ioverlap = 0 in Record 1 is required.
c
c=====Record 3: wfnfile
c
c        Name of the .WFN file containing the definition of the wave
c        function. In the case of a CASSCF calculation this includes 
c        the determinant coefficients and the occupation numbers of 
c        their molecular spin-orbitals. Files of this type are presently
c        available for .WFN files coming from a domestic version of the 
c        GAMESS10 code. It is clear that the overlap integrals read 
c        in from 'filedat' must correspond to the wave function given 
c        in 'wfnfile'
c
c-----After record 3 the following input lines are allowed.
c     (records 4 and 5.i are mandatory and must be consecutive)
c
c=====Record 4: NGROUP ngroup 
c
c        The number domains in which the molecule is divided.
c
c        If ngroup is less than or equal to 0, or an integer is not given
c        after the NGROUP keyword, each basin is taken as an independent 
c        domain with a minimum (maximum) electronic population equal to 0 
c        (N_A) for ALPHA electrons and 0 (N_B) for BETA electrons, where 
c        N_A and N_B are the total number of ALPHA and BETA electrons in 
c        the molecule. In this case Record 5.i must be skipped.
c
c        The same applies if the 'NGROUP' keyword is not present in the 
c        input file.
c
c        If 'ngroup' is given a negative value the 'i' variable in the
c        following record runs for 1...NGROUP-1, and last group is 
c        intended to be formed by all those atoms not included in the
c        first NGROUP-1 groups.
c
c=====Record 4bis: NFRAG nfrag
c
c        This is equivalent to the order 'NGROUP ngroup' with the 
c        difference that 'nfrag' can be much larger than 'ngroup'
c        because the 'nfrag' value is exclusively used in the 
c        'nfnoeos.f' and 'lewisoqs.f' routines that compute the 
c        effective oxidation state (EOS) of the defined fragments
c        and the Lewis structures, respectively. If this order is 
c        not given in the input file of EDF but the 'nfnoeos.f' or
c        'lewisoqs.f' routine is invoked with the 'oqseos' and
c        'lewisoqs' orders, respectively (see below), 'nfrag' is
c        taken as equal to 'ngroup'.
c
c=====Record 5.i (i=1,NGROUP)
c
c         Record 5.i may have two different formats, (1) and (2). 
c
c         Format (1) is:
c
c    (1)  nfugrp(i), [ifugrp(j,i),j=1,nfugrp(i)], 
c              minelecA, maxelecA, minelecB, maxelecB
c
c         where
c
c         nfugrp(i)     = Number of atoms in group i.
c         ifugrp(j,i)   = Indices of atoms that belong to group i.
c         minelecA      = Minimum number of ALPHA electrons in group i.
c         maxelecA      = Maximum number of ALPHA electrons in group i.
c         minelecB      = Minimum number of BETA  electrons in group i.
c         maxelecB      = Maximum number of BETA  electrons in group i.
c
c         Format (2) is:
c
c    (2)  nfugrp(i),[ifugrp(j,i),j=1,nfugrp(i)]
c
c         In this case, EDF associates to group 'i' the values
c
c         minelecA = 0   maxelecA = NALPHA
c         minelecB = 0   maxelecB = NBETA
c
c         where NALPHA and NBETA are the total number of ALPHA and BETA
c         electrons of the molecule, as computed in the rdwfn.f routine.
c
c         In both cases, minelecA, maxelecA, minelecB, maxelecB MUST 
c         INCLUDE the CORE electrons, possibly defined in the COREMO
c         order below, i.e. minelecA, maxelecA, minelecB, maxelecB DO 
c         NOT refer exclusively to valence electrons (not equal to the 
c         total number of electrons when mocore # 0), but to total 
c         electrons. EDF assumes that each core MO accounts of one ALPHA 
c         electron and one BETA electron.
c
c=====Record 5bis.i (i=1,NFRAG)
c
c         These 5bis.i orders are equal to the above 5.i with the 
c         difference that 'i' in this case goes from 1 to NFRAG, whereas 
c         'i' in the 5.i orders goes from 1 to NGROUP. Its format is
c
c         nwithin(i),[ifrag(j,i),j=1,nwithin(i)]
c
c         where ifrag(1:nwithin(i),i) are the indices of the atoms
c         that belong to the fragment 'i'.
c
c=====Record 6: AOMNORM atom
c
c         When the keyword AOMNORM is read in, the program modifies
c         the elements of the atomic overlap matrix corresponding
c         to the center 'atom' such that SUM (i=1,ncent) aom(i,m,j) 
c         will be exactly 0.0 for m.ne.j or 1.0 for m.eq.j. If 'atom'
c         is not read in the program takes 'atom'=last atom of the 
c         input .WFN file.
c
c         After the normalizacion process is carried out, the new AOM 
c         elements are written in the file filedat//'.normalized'
c
c
c=====Record 7: EPSDET real-number 
c
c        The real-number has the following meaning: A product Ci * Cj
c        of a pair of determinants i and j in a CASSCF calculation is 
c        neglected when abs(Ci*Cj) < real-number. Default is 0.0.
c
c=====Record 8: EPSPROBA real-number 
c
c        The real-number has the following meaning: The computation 
c        of any probability for BETA electrons is avoided if the 
c        highest probability for ALPHA electrons is smaller than
c        this 'real-number'
c
c=====Record 9: PROBCUT real-number  
c
c        In the output of probabilities, only those probabilities 
c        greater than PROBCUT are written. Default is 0.0
c
c=====Record 10: TOLAOM real-number  
c
c        If the sum over all the centers of an atomic overlap matrix
c        element, SUM_n AOM(n,i,j), differs from delta_ij by more than
c        the 'real-number' the program stops since this means that
c        the AOM has not been appropriately computed. The default value
c        of real-number is 0.01.
c
c=====Record 11: NDETS ndets  
c 
c        In a CASSCF calculations only the first ndets determinants in 
c        the wavefunction are used (the wavefunction is renormalized).
c        If this record is not given EDF uses all the determinants of
c        the wave function.
c
c=====Record 12: EPSWFN epswfn
c
c        This is another order, related to the previous one, that can be
c        used to control the number of determinants of the input .WFN 
c        file that is actually used to expand the wavefunction. All the 
c        determinants whose (absolute value of the) coefficient in the 
c        wavefunction is lesser than 'epswfn' are neglected. After this, 
c        the wave function is renormalized. The default value is 0.0.
c
c=====Record 13: EPSWFNLOC epswfnloc
c
c        Equivalent to the above order but with the WFN expanded in
c        terms of (isopycnic) localized MOs. This order is relevant 
c        only if the FULLOC order is active (see below).
c
c=====Record 14: LARGE or VERBOSE or BIGOUTPUT
c
c        Large output is requested. The default is short output.
c
c=====Record 15: RECUR
c
c        By default, when the number of groups is 2, the recurrence
c        relations of Cançes et al are not used. With this order, 
c        such recurrences are employed even in the case of MD WFN's.
c
c
c=====Record 16: FULLOC [ ALLCENTER ]
c
c        Localize the input natural orbitals. By default the isopycnic
c        localization is not carried out. If the 'allcenter' keyword is
c        not given the localization is carried out using the groups
c        defined in the input to contruct the localization function 
c        which is maximized. On the contrary, if 'allcenter' keyword is
c        is given the localization function is built using every center 
c        as an individual fragment. 
c
c=====Record 17: NOPROB
c
c        If this keyword is given the exact EDF using the wavefunction 
c        defined in the input .WFN file is not computed. The defaut is 
c        to compute it.
c
c=====Record 18: PROBLOC
c
c        If this keyword is given, the EDF using the wavefunction
c        expressed in terms of the Localized MOs is computed. By default
c        this EDF is not computed. Clearly, this keyword is relevant only
c        in case that a previous localization of the input natural MOs 
c        has been performed. Since isopycnic localized MOs are normalized
c        but not orthogonal, the sum of all probabilities of the EDF in 
c        this case is not, in general, equal to 1.0.
c
c=====Record 19: NOPROBX
c
c        If this keyword is given, the approximate EDF using the wave-
c        function expressed in terms of the Localized MOs is not computed.
c        The defaut is to compute it. This keyword is relevant only
c        in case that a previous localization of the input natural MOs 
c        has been performed. The approximation in computing the EDF 
c        consists in assuming that all the localized MOs are orthogonal 
c        within the different fragments in which the molecule has been 
c        divided.
c
c=====Record 20: COREMO overlap-canonical
c
c        This keyword indicates that any molecular orbital of the input
c        wavefunction with a self-overlap in any of the fragments greater
c        than 'overlap-canonical' is considered as a fully-localized 
c        (i.e. core) orbital in that fragment. Consequently, it is 
c        excluded from the computation of the EDF. If the keyword 
c        COREMO is given but the 'overlap-canonical' value is not read,
c        the program assumes overlap-canonical = 0.99.
c
c=====Record 21: NOORDER
c        This keyword indicates that spinless probabilities are not 
c        ordered by decreasing value in the output file. By default, 
c        they are ordered provided that its number is smaller that the
c        NSTACK parameter, currently set to 4000.
c
c=====Record 22: DELETE n_1 n_2 n_3 ...
c        Let us assume that a multideterminantal wavefunction with NCORE
c        MOs has been read in. With this order, provided n_i <= NOCORE, 
c        the ith MO is elliminated from the wavefunction. After this
c        order, the actual number of molecular orbitals and CORE orbitals
c        (both ALPHA and BETA) are decreased by the number of valid n_i's 
c        which are read in this order.
c
c=====Record 23: SUPPRESS n_1 n_2 n_3 ...
c        Suppress from the AOM file the MOs n_1, n_2, n_3,..., writing
c        a new AOM file with the AOM integrals between the non-deleted
c        MOs, and the MOs n_1, n_2, n_3 from the WFN file, writing a 
c        new WFN file with the rest of MOs. The last task is done
c        in case of single- or multi-Det wavefunctions, although in the
c        last case de routine is not fully tested yet. After performing 
c        these tasks the program stops.
c
c=====Record 23: LOCORE
c        Localize the CORE MOs from a multideterminantal wave function,
c        and replaces the canonical MOs of the wavefunction by these 
c        localized MOs. From this point on the original MOs are lost.
c
c=====Record 24: LOCSOME n1,n2, ...
c        Performs an isopycnic localization on the subset of canonical
c        MOs given by n1,n2, .... These can be typically CORE orbitals
c        in a multideterminantal wavefunction. After the localization
c        the AOM is appropriately modified.
c
c=====Record 25: NOMEM
c        For closed-shell systems and multiDet WFNs the ALPHA and BETA 
c        probabilities are saved in arrays instead of being stored in 
c        temporary files. This is the default (memedf=.true.) and in 
c        this case the routine 'mcalcedf.f' instead of 'calcedf.f' is 
c        used. This keyword is used to change this option (i.e. to 
c        set memedf=.false.)
c
c=====Record 26: PRSRS N_1 N_2 ... N_G
c        Compute the probability of the real space resonance structure
c        (PRSRS) in which N_1 electrons are in group 1, N_2 electrons in
c        group 2, ..., and N_G electrons in the last group. It is 
c        possible to put as many lines as this one in the input file.
c
c=====Record 27: MAXPOP maxgroup_1 maxgroup_2  ...
c
c        Define the maximum number of ALPHA or BETA electrons in each
c        group when the PRSRS oder is applied. Only the spin-resolved
c        resonance structures for which the number of ALPHA or BETA 
c        electrons in group "i" is smaller than or equal to maxgroup_i 
c        are considered. This order must be given after the PRSRS order.
c
c=====Record 28: MINPOP mingroup_1 mingroup_2  ...
c
c        Define the minimum number of ALPHA or BETA electrons in each
c        group when the PRSRS oder is applied. Only the spin-resolved
c        resonance structures for which the number of ALPHA or BETA 
c        electrons in group "i" is greater than or equal to mingroup_i 
c        are considered. This order must be given after the PRSRS order.
c
c=====Record 29: RANDOM minrandom maxrandom
c
c        Define the mininum (minrandom) and maximum (maxrandom) values
c        of a real*8 random number to be used in solving the linear
c        system of equations to obtain the EDF probabilities in the
c        'binedf.f' routine.
c
c=====Record 28: TWOCENDI [ epsbond ]
c
c        Activates the computation of all the two-center delocalization
c        indices. 'epsbond' is a threshold value such that atoms A and B
c        are considered to be bonded if DI_AB > epsbond. The default va-
c        lue of epsbond is 1D0. This order can only be used with Single
c        Determinant Wavefunctions (SDW). Ten calculations are perfor-
c        med usig epsbond, epsbond/2, ..., epsbond/2**9.
c
c=====Record 29: KPOWER pow
c
c        K power of the iterative Becke integration scheme. Default
c        is pow=3.
c
c=====Record 30: LEBEDEV nang
c
c        Number of angular points in Lebedev integrations. Default is
c        nang=434
c
c=====Record 31: NRAD nrad
c
c        Number of radial points in Becke integration. Default is 
c        nrad=100
c
c=====Record 32: IRMESH irmesh
c
c        Mapping of the r coordinate in Becke integration. 
c        Default is irmesh=1
c
c=====Record 33: RHOPOW rhopow
c
c        Power of each rho_A when computing w_A as 
c        w_A=rho_A^rhopow/SUM(rho_A^rhopow) in the PROMRHO partition.
c
c=====Record 34: CANLOC [ ALLCENTER ]
c
c        Localize the input canonical orbitals. By default the isopycnic
c        localization is not carried out. If the 'allcenter' keyword is
c        not given the localization is carried out using the groups
c        defined in the input to contruct the localization function 
c        which is maximized. On the contrary, if 'allcenter' keyword is
c        is given the localization function is built using every center 
c        as an individual fragment.
c
c=====Record 35: CHEMLOC
c
c        Idem as the above one but using the so-called chemical loc.
c
c=====Record 36: AOMLOC
c
c        Idem as the above one but using another variant of the 
c        chemical localization.
c
c=====Record 37: DAFHLOC
c
c        Idem as the above one but using the diagonalization of the
c        DAFH instead of the AOM.
c
c=====Record 38: OPENLOC
c
c        Idem as the above one but using the diagonalization of the
c        Open Quantum Sysmtem (OQS) 1-RDMs of atoms, pairs of atoms,
c        etc, i.e S_a^{1/2} * 1-RDM * S_a^{1/2}.
c
c=====Record 39: DAFHSELECT n n_1 n_2 ... n_n
c
c        When this record is read in, the search of localized nc-2e MOs
c        IS NOT carried out for all the n-tuples with a level 'n' but
c        only for the n-tuple of atoms indicated by 'n_1 n_2 ... n_n'. 
c        There can be in the input file as many 'DAFHSELECT' orders as 
c        desired. If there are no DAFHSELECT orders with a given value 
c        of 'n', then ALL n-tuples of atoms with 'n' atoms are explored. 
c        On the other hand, if a given value of 'n' is skipped in the order 
c        'MXCENT' (see below), then all DAFHSELECT orders with this value
c        of 'n' are irrelevant. Let us consider an example. Assume a
c        molecule with 5 atoms, and consider that the following lines
c        are present in the input file:
c
c        DAFHSELECT 3 1 2 3
c        DAFHSELECT 3 1 3 5
c        DAFHSELECT 4 1 2 3 4
c
c        Then, the following nc-2e searches for localized MOs are performed:
c
c        All 1c-2e searches
c        All 2c-2e searches
c            3c-2e searches for the 3-tuples '1 2 3' and '1 3 5'
c            4c-2e searches for the 4-tuples '1 2 3 4' 
c        ALl 5c-2e
c        ALl 6c-2e
c           "
c           "
c        The above list goes up to a value of 'n' equal to 'mxcent' 
c        (see the MXCENT order below) and the values of 'n' equal to 
c        'n1', 'n2', ..., in this order are skipped.
c
c        At this moment, this restricted searchs works only in the case
c        of localized nc-2e MOs of Single-Det WFNs within the DAFHLOC
c        analysis.
c
c=====Recors 40: DODAFH
c                 :
c                 :
c                 :
c                END
c
c        The aim of this order is similar to that of the DAFHLOC. The 
c        difference is that DODAFH has a full control of the DAFH 
c        localizations actually performed. Within the DODAFH ... END
c        environment one puts as many orders like this one:
c
c        ncen cutoff (icen(j),j=1,ncen) [DEPLET]
c
c        where 'ncen'   ---> Indicates a localization of NCEN centers
c              'cutoff' ---> The value (close to but smaller than 1.0)
c                            to discriminate whether a found MO is 
c                            localized or not.
c              'icen()' ---> The indices of the NCEN atoms.
c              DEPLET   ---> means that the full DAFH is depleted after
c                            the last tried localization.  
c
c        The following is an example for the C_2H_6 molecule
c
c        dodafh
c          1 0.95 1
c          1 0.95 2   deplet
c          2 0.95 1 2 deplet
c          2 0.95 1 5
c          2 0.95 1 6
c          2 0.95 1 7 deplet
c          2 0.95 2 3
c          2 0.95 2 4
c          2 0.95 2 8 deplet
c        end
c
c=====Recors 41: OTHEROPEN
c                 :
c                 :
c                 :
c                END
c
c        Idem as the above one but using the diagonalization of the
c        Open Quantum Sysmtem (OQS) 1-RDMs of atoms, pairs of atoms,
c        etc, i.e S_a^{1/2} * 1-RDM * S_a^{1/2}.
c
c=====Record 40: CRITICS critics
c
c        When CANLOC or CHEMLOC is invoked, the routine determines the set
c        of localized MOs that are partially localized in each atom or 
c        group. A MO fulfils this condition is S_ii^A > critics. The
c        default value is critics=0.02
c
c======Recors 41: COVX covx
c
c     If the distance between two atoms A and B is smaller than
c     (RA+RB)*COVX, a bond is assumed to exist between A and B.
c     RA and RB are the covalent radius stored in the array covrad()
c     of the 'connect.f' routine. The default value is 1.2.
c
c
c======Recors 42: DAMPS damps
c     
c     This order is relevant only for the 'CHEMLOC' localization.
c     In every cycle of localization the value of each CRITOV(1..mxcent)
c     is multiplied by DAMPS (a number that should be smaller than 1.0) 
c     with the aim of reducing the localization requirement. 
c     The default value is 0.95
c
c
c======Record 43: MXCENT mxcent [ n1 n2 ... ]
c
c     In the progressive localization of MOs based on the diagonalization
c     of the DAFH in 1-, 2-, 3-,... centers, 'mxcent' defines the maximum
c     number of centers actually explored. The default value is 6. If the
c     cases with 'n1', 'n2', ... centers are to be skipped, they should
c     be indicated after the value of 'mxcent'. For instance, the order
c    
c     MXCENT 6 3 4 
c
c     analyzes the 1- 2- 5- and 6-center localizations.
c
c======Record 44: CRITOV [ { nn | critval } ]
c
c     If the overlap of an MO with itself in a n-center domain is greater 
c     than  critov(n), the routines 'chemloc', 'itaomloc', and 'itdafh' 
c     consider that this MO is already localized on that domain. 
c     The default values for all CRITOV(1..mxcent) is 0.95. If no 
c     parameters are given after the CRITOV keyword none of the current
c     values of CRITOV(1..mxcent) are changed. If an integer number 'nn'
c     such that (1.LE.NN.LE.MXCENT) is given after the CRITOV keyword
c     the new CRITOV(nn) is given the value 'critval' (in case a real
c     number is read in after reading the value of 'nn'). Finally,
c     if a real number 'critval' is read in after the CRITOV keyword,
c     all CRITOV(1..mxcent) values will have the value 'critval' from 
c     now on. When the 'DAFHLOC' order is invoked and the WFN file
c     corresponds to a correlated wavefunction the critiv() values
c     are not checked againts the diagonal AOM values but againts the
c     eigenvalues of the pertinent DAFH matrix.
c
c======Record 44a: CRITOV1  critval 
c
c      Sets critov(1) equal to 'critval'. This order is equivalent
c      to CRITOV 1 critval
c
c======Record 44b: CRITOV2  critval 
c
c      Sets critov(2) equal to 'critval'. This order is equivalent
c      to CRITOV 2 critval
c
c=====Record 45: MAXBOND maxbond
c
c     In the progessive DAFH density partitioning (when pairs of atoms
c     are being analized) 'maxbond' represents the maximum number of 
c     bonds necessary to travel from 'atom1' to 'atom2' which is actually
c     explored; i.e. if going from 'atom1' to 'atom2' requires more
c     tham 'maxbond' bonds the pair 'atom1-atom2' is not explored.
c     The default value of 'maxbond' is 1, i.e. only directly bonded
c     atoms are considered. The maximum value of 'maxbond' is 
c     'maxmaxb = 4'.
c
c=====Record 47: SKIPH
c
c     In the progessive DAFH density partitioning (when single atoms
c     are being analized) when 'skiph=.true.' Hydrogen atoms are skipped.
c     Otherwise, hydrogen atoms are explicitely taken into account.
c     The default is 'skiph=.false.'.
c
c=====Record 48: EPSNEG [epsneg]
c
c     In the OPENLOC or OTHEROPEN localizations, diagonalizations
c     are involved of definite positives matrices that should have
c     all of their eigvenvalues positive. Numerical errors, however,
c     make that this not necessarily happens and some of the 
c     eigenvalues can be marginally negative. The value of 'epsneg'
c     is the largest value (in absolute value) allowed in order
c     that the method ends nicely. For instance, if epsneg=-1D-6
c     and all the eigenvalues are greater that this number the 
c     process continues up to the end. However, if an eigvenvalue
c     is smaller that -1D-6, for instance -1D-5, the program aborts.
c     The default value of epsneg is -1D-6. This value can be
c     changed with the present order.
c
c=====Record 49: EPSEIGEN [epseigen]
c
c     In the diagonalization of [S^A]^(1/2) * 1RDM * [S^A]^(1/2), the
c     eigenvalues smaller than 'epseigen', as well as the associated
c     eigenvectors, correspond to FNOs mostly localized outside the 
c     fragment A and will ge ignored after the above diagonalizaition 
c     is performed.
c
c=====Record 50: SMALLEIGEN [smalleigen]
c
c-----In the computation of [S^A]^(-1/2) a division by the square root
c     of an eigenvalue of S^A is involved. In many cases, such eigenvalue
c     is zero or very small. To avoid a division by 0.0 the eigenvalue
c     of S^A is replaced by 'smalleigen' if it is smaller than this number.
c
c=====Record 51: OFMO
c
c        Orthogonalize exactly the canonical and natural MOs using
c        the OFMO method.
c
c=====Record 52: DOSYMMETRY [ { 1 | 2 | 3 } ]
c
c        Performs symmetry analysis of the molecule with method 1,
c        method 2, or method 3. The default is method 1.
c
c        1 ---> Symmetry analysis with the standard code of promolden
c        2 ---> Symmetry analysis with O. Beruski and L.N. Vidal code 
c               ( lnvidal@utfpr.edu.br )
c        3 ---> Symmetry analysis with Jose Luis Casals Sainz code
c
c=====Record 53: TOLSYM [tolsym]
c
c        With this order we change the value of the 'tolsym' parameter
c        that is used in symmetry modules and measures the tolerance
c        for identity of points,matrices,... Default value is 1e-6.

c=====Record 54: EPSSYM [epssym]
c
c        With this order we change the value of the 'epssym' parameter
c        that determines if molecular orbitals in 'transor.f' routine
c        are considered linearly independent or not. Default value is
c        1e-6.

c=====Record 55: EPSORB [epsorb]
c
c        With this order we change the value of the 'epsorb] parameter
c        that determines if molecular orbitals in 'transor.f' routine
c        break symmetry or not. Default value is 1e-4.
c
c=====Record 56: POINTGROUP pntgrp
c
c        With this record we specify the point group of the molecule
c        instead of determining it within the symmetry routines.
c
c=====Record 57: NOSYMMETRY
c
c        The symmetry is not used in the 'sym.f' module, i.e. the
c        point group of the molecule is assumed to be C1
c
c=====Record 58: OPENSYS
c
c        With this order, we perform an open-quantum-system (OQS) analysis, 
c        with the aim of determining the Oxidation State (OS) of an atom
c        in the molecule, in the spirit of the EOS method of Salvador et
c        al. Before this order, a previous calculation of the alpha and
c        beta 1-RDMs is necessary. In case of a CCWFN, only closed-shell
c        cases are possible to analyze. After the OPENSYS order has been
c        completed, EDF stops. 
c AQUI
c        determining
c        the DNOs of all the GROUPS in a way, alternative to the DAFH 
c        analysis. This is done only in case that the wavefunction is not
c        of CCIQA type.
c
c=====Record 60: DORDM
c
c        Compute the 1-RDM of the molecule in case its wavefunction is not 
c        of a CCIQA type.
c
c=====Record 61: BONDING
c                   :
c                   :
c                ENDBONDING
c
c         The 'BONDING' keyword must be put after the 'NGROUP' order
c         defining the fragments of the molecule.
c
c         With the set of orders between the BONDING and ENDBONDING 
c         keywords, we try to represent the full EDF of the molecule 
c         (partitioned into a given set of fragments) as a direct product 
c         of [2c,2e] [and possibly (3c,2e)] chemical links. After
c         the 'BONDING' keyword the following orders can be given that
c         are exclusively relevant of the 'optedf.f' routine:
c
c        1) ENDBONDING
c        2) TYPE
c        3) TYPE3C
c        4) PAIR
c        5) TRIO
c        6) EPSBOND
c        7) INICVAR
c        8) PRIN
c
c        9) PLUS/MINUS
c
c       10) WEDF
c       11) WPOP
c       12) WDIS
c
c     The meaning of all of these orders is explained below.
c
c     1) 'ENDBONDING'
c
c     Sintax: ENDBONDING
c
c     Meaning: This order means that all the input that the routine 
c     'optedf' needs to work has been given, so that the calculation 
c     can start. This keyword in the input file is MANDATORY.
c
c     2) 'TYPE'
c
c     Sintax: TYPE type qq(type) ifixq(type) ffty(type) ifixf(type)
c
c     Meaning: This order defines a type of (2c,2e) bond. The maximum
c     value of 'type' is MAXBOND, a parameter whose current value is
c     maxbond=18. To increase the maximum possible number of types of
c     (2c,2e) bonds just increase the MVAR, MVARH, and MAXBOND values 
c     in the sentence 'parameter (mvar=36,mvarh=mvar/2,maxbond=18)' 
c     below and recompile.
c
c     'qq(type)' is 'the charge value of the bond', which is a measure
c     of the polarity of the bond. Its value is optimized in case that
c      the read in value of 'ifixq(type)' is non zero. Otherwise, the
c     value of 'qq(type)' is fixed. 
c
c     'ffty(type)' if the 'correlation factor of the bond'. Its value 
c     is optimized in case that the read in value of 'ifixf(type)' is
c     non zero. Otherwise, the value of 'ffty(type)' is fixed.
c
c     NO default values exist for 'qq(type)' and 'ffty(type)' and the
c     five values after the 'TYPE' keyword are mandatory. Tipical input
c     values for 'qq(type)' and 'ffty(type)' are 0.5 and 0.0,
c     respectively.
c
c     3) 'TYPE3C'
c
c     Sintax: TYPE3C type p200 p020 p002 p110 p101 f200 f020 f002 f110 f101
c
c     Meaning: This order defines a type of (3c,2e) bond. The maximum
c     value of 'type' is MAXBOND3C, a parameter whose current value is
c     maxbond3c=6. To increase the maximum possible number of types of
c     (3c,2e) bonds just increase the MAXBOND3C and MVAR3C values in the 
c     sentence 'parameter (mvar3c=30,maxbond3c=6)' below and recompile. 
c
c     'p200 p020 p002' are the input values for the probability that 
c     both electrons of the 3c bond are located in the first center (p200),
c     in the second center (p020), or in the third center (p002).
c
c     'p110' is the input value for the probability that one of
c     the electrons of the 3c bond is located in the first center and 
c     the second electron in the second center.
c
c     'p101' is the input value for the probability that one of
c     the electrons of the 3c bond is located in the first center and 
c     the second electron in the third center.
c
c     'f200 f020 f002 f110 f101' are five integer numbers that determine
c     if 'p200 p020 p002 p110 p101' are fixed or optimizable parameters.
c     For instance, f200=1 or 0 means that 'p200' is fixed or optimizable
c     respectively.
c
c     COMMENT: Fixing parameters does not work very well yet.
c
c      4) 'PAIR'
c
c      Sintax: PAIR fragment1 fragment2 type [ rtype ]
c
c      Meaning: This order indicates that a (2c,2e) bond of type 'type'
c      purportedly exists between the groups of atoms 'fragment1' and 
c      'fragment2'. 

c      'rtype' indicates the resonance structure used. The default is
c      rtype=1. It is recommended at this moment NOT TO ENTER THE
c      rtype VALUE IN THIS ORDER.
c
c     5) 'TRIO'
c
c     Sintax: TRIO fragment1 fragment2 fragment3 type [ rtype ]
c
c     Meaing: This order indicates that a (3c,2e) bond of type 'type'
c     purportedly exists between the groups of atoms 'fragment1',
c     'fragment2', and 'fragment3'.

c     'rtype' indicates the resonance structure used. The default is
c     rtype=1. It is recommended at this moment NOT TO ENTER THE
c     rtype VALUE IN THIS ORDER.
c
c     6) 'EPSBOND'
c
c     Sintax: EPSBOND epsbond
c
c     Convergence threshold in the EDF fitting. This order is not nece-
c     ssary, since a default value for 'epsbond' exists [10^(-4)].
c     
c
c     7) 'INICVAR'
c
c     Sintax: INICVAR h0
c
c     Meaning: 'h0' is a parameter used in the EDF fitting. Not necessary 
c     since a default value exists (0.1).
c
c     8) PRIN prin
c
c     Sintax: PRIN prin
c
c     Meaning: 'prin' is a printing level in the EDF fitting. Not necessary
c     since a default value exists (0).
c
c     9) PLUS/MINUS
c
c     Sintax { PLUS | MINUS } group  Resonance_structure
c
c     Meaning:
c
c
c
c
c
c     10) WEDF
c
c     Sintax: WEDF [ weigedf ]
c
c     Meaning: The N-electron EDF in terms of (2c,2e) [and possibly (3c,2e)]
c     links is optimized by minimizing:
c
c     delta = WEDF * A + WPOP * B + WDIS * C
c
c     where WEDF, WPOP, and WDIS must be real positive numbers,
c
c     A = SUM_{i} [p^approx(i)-p^exact(i)]^2,
c
c     where 'i' runs over all computed real space resonance structures (RSRS),
c
c     B = SUM_{i} [N^approx(i) - N^exact(i)]^2,
c
c     where 'i' runs over all the fragments in which the molecule has been
c     divided and N^approx(i) and N^exact(i) are the approximated and exact
c     electron population of group 'i',
c
c     C = SUM_{i>j} [ delta^approx(i,j) - delta^exact(i,j}]^2,
c
c     where i>j runs over all pairs of groups and delta^approx(i,j) and
c     delta^exact(i,j} are the approximate and exact delocalization
c     indices between the groups 'i' and 'j'.
c
c     The default values of WEDF, WPOP, and WDIS are 1.0, 0.0, and 0.0,
c     respectively. With the 'WEDF [ weigedf ]' order a new value is
c     assigned to the weight of quantity 'A'
c
c     11) WPOP 
c
c     Sintax: WPOP [ weigpop ]
c
c     Meaning: A new value is assigned to the weight of quantity 'B'
c
c     12) WDIS
c
c     Sintax: WDIS [ weigdis ]
c
c     Meaning: A new value is assigned to the weight of quantity 'C'
c
c-----------------------------------------------------------------------
c
c=====Record 62: ROHFLOC
c
c     Localizes separately the ALPHA and BETA spin-orbitals of a ROHF WFN 
c     and stops. A WFN file with the localized MOs is written. An AOM file
c     with the overlaps between the localized MOs is also written.
c
c=====Record 63: UHFLOC
c
c     The same as in the previous ROHFLOC order but with an UHF WFN.
c
c=====Record 64: DOENTROPY
c
c     Mutual Entropy Information (MEI) is obtained after the EDF is 
c     computed. By default the MEI is not obtained since it usually 
c     implies a very large output.
c
c=====Record 65: DIRPROB
c
c     This order computes by a pure brute force method the probabilities
c     of the sub-set of RSRSs with electron populations between minpopul(i) 
c     and maxpopul(i) values. At this time, it only works for SDWs. IT IS 
c     HIGHLY RECOMMENDED NOT TO USE this order in molecules with many
c     electrons and/or divided into a large number of fragments.
c
c=====Record 66: FNO atom1 atom2 ...
c
c     Determines the Fragment Natural Orbitals of the fragment of the 
c     current molecule formed by atoms 'atom1', 'atom2', ... There can 
c     be as many FNO orders as desired in the input file. After completing
c     this order, EDF stops.
c
c=====Record xx: CGEDF
c=====Record xx: RHOCOND
c=====Record xx: ALLRHOC
c=====Record xx: RHOALLS
c=====Record xx: NINDEP ngroup
c=====Record xx: MAXWFN
c
c=====Record 6x, first ascii character = '#'
c
c        This indicates a non-executable comment line.
c
c=====Record 6x: END
c
c        This indicates the end of the input of the current calculation
c 
c.......................................................................
c
      USE          space_for_wfnbasis
      USE          space_for_wfncoef
      USE          space_for_cidet
      USE          mod_sym
      USE          mod_periodic
      USE          space_for_sym
      USE          space_for_rdm1
      USE          space_for_sgarea
      USE          space_for_aomspin
      include     'implicit.inc'
      include     'param.inc'
      include     'wfn.inc'
      include     'corr.inc'
      include     'stderr.inc'
      include     'constants.inc'
      include     'lengrec.inc'
      include     'fact.inc'
      include     'sym.inc'
      include     'datatm.inc'
      include     'integ.inc'
      include     'mline.inc'
      integer(kind=4), parameter ::   numrsrs =   50
      integer(kind=4), parameter ::   maxexp  =  100
      integer(kind=4), parameter ::   lline   = 2000
      real(kind=8),    allocatable,dimension (:)     :: eos,eosab
      real(kind=8),    allocatable,dimension (:,:)   :: vrot
      real(kind=8),    allocatable,dimension (:)     :: work
      real(kind=8),    allocatable,dimension (:,:,:) :: aom,aomx
      real(kind=8),    allocatable,dimension (:,:,:) :: aomnat
      real(kind=8),    allocatable,dimension (:,:,:) :: aomtopmod
      real(kind=8),    allocatable,dimension (:,:)   :: cc
      real(kind=8),    allocatable,dimension (:,:)   :: v1mata
      real(kind=8),    allocatable,dimension (:,:)   :: v1matat
      real(kind=8),    allocatable,dimension (:)     :: d1mata
      real(kind=8),    allocatable,dimension (:,:,:) :: dafh
      real(kind=8),    allocatable,dimension (:,:,:) :: sgallc
      real(kind=8),    allocatable,dimension (:,:)   :: canmo
      real(kind=8),    allocatable,dimension (:,:)   :: sgtmp
      real(kind=8),    allocatable,dimension (:)     :: vtmp
      real(kind=8),    allocatable,dimension (:)     :: oc
      integer(kind=4), allocatable,dimension (:,:)   :: ifugrp
      integer(kind=4), allocatable,dimension (:,:)   :: ifrag
      integer(kind=4), allocatable,dimension (:,:)   :: eleca,elecb
      integer(kind=4), allocatable,dimension (:)     :: nfugrp
      integer(kind=4), allocatable,dimension (:)     :: nwithin
      integer(kind=4), allocatable,dimension (:)     :: iaom
      integer(kind=4), allocatable,dimension (:)     :: ipvt
      integer(kind=4), allocatable,dimension (:)     :: icore,ival
      integer(kind=4), allocatable,dimension (:)     :: icorea,icoreb
      integer(kind=4), allocatable,dimension (:)     :: ivaloc
      integer(kind=4), allocatable,dimension (:,:)   :: icogrp
      integer(kind=4), allocatable,dimension (:)     :: mocogrp
      integer(kind=4), allocatable,dimension (:)     :: mocogrpa
      integer(kind=4), allocatable,dimension (:)     :: mocogrpb
      integer(kind=4), allocatable,dimension (:)     :: isup,isupx
      integer(kind=4), allocatable,dimension (:)     :: ordsup
      integer(kind=4), allocatable,dimension (:)     :: maxpopul
      integer(kind=4), allocatable,dimension (:)     :: minpopul
      integer(kind=4), allocatable,dimension (:)     :: maxalpha
      integer(kind=4), allocatable,dimension (:)     :: minalpha
      integer(kind=4), allocatable,dimension (:)     :: maxbeta
      integer(kind=4), allocatable,dimension (:)     :: minbeta
      integer(kind=4), allocatable,dimension (:,:)   :: resnca,resncb
      integer(kind=4), allocatable,dimension (:,:)   :: inside
      integer(kind=4), allocatable,dimension (:)     :: ninside
      real(kind=8),    allocatable,dimension (:,:)   :: eta
      integer(kind=4), allocatable,dimension (:)     :: atnum
      real(kind=8)     rcond,deter(2)
      real(kind=8)     aom1,aom2,aomtot,tmpsg
      integer(kind=4), parameter :: maxcritov = 10
      real(kind=8),    parameter :: critovdef = 0.95D0
      real(kind=8)     critov(maxcritov)
      logical          ixcent(maxcritov)
      logical          dafhselect(maxcritov)
      integer(kind=4)  ndafhsel   (maxcritov)
      integer(kind=4)  idafhsel   (maxcritov,maxcritov,maxcritov)
      integer(kind=4)  stdout,stdin,udat,udatnw,leng,lw
      logical          setint,setword,setdble,memedf
      logical          ok,ok1,ok2,ok3,ok4,ok5,ok6
      logical          mulliken,lowdin,qtaim,mindef,becke,allcloc
      logical          mindefrho,netrho,promrho,hesel,notaom,skipcan
      logical          skiploc,skiplocx,largwr,short,wfncut,okaom,goon
      logical          iscore,ditwocent,nindep,ngless,localize
      logical          coreloc,corecan,orderp,locore
      logical          canloc,chemiloc,dafhloc,aomloc
      logical          bonding,maxwfn,aomnorm,recur,dosymmetry 
      logical          warn,warno,skiph,exfil,otherdafh
      logical          opensys,ok1rdm,oqseos,dolewis
      logical          openloc, otheropen, dofno, doentropy, rohfloc
      logical          uhfloc, dodirprob, noloc, lewis, defragoqs
      logical          salvanu,scuseria
c
      character*(mline) filedat,cicoef,word,line,uppcase,wfnfile
      character*(mline) fileinp
      character*(lline) largeline,oldlarge
      character*(mline) fileloc,nline
      character*(mline) wfnloc,typepart,nres(numrsrs)
      character*(mline) lower
      character*16     namebas(100)
      integer(kind=4)  pow,nrad,nang,irmesh
      real(kind=8)     ovcritloc,ovcritcan
      real(kind=8),    parameter ::   highover      =  0.99d0
      real(kind=8),    parameter ::   smallover     =  0.80d0
      real(kind=8),    parameter ::   defepswfn     =  0.00d0
      real(kind=8),    parameter ::   deftolaom     =  0.05D0
      real(kind=8),    parameter ::   defcritics    =  0.02D0
      real(kind=8),    parameter ::   defahesel     =  2.00D0
      real(kind=8),    parameter ::   defbhesel     =  2.00d0
      real(kind=8),    parameter ::   defdamps      =  0.95d0
      real(kind=8),    parameter ::   defcovx       =  1.20d0
      real(kind=8),    parameter ::   ranmindef     = -5.00d0
      real(kind=8),    parameter ::   ranmaxdef     = +5.00d0
      real(kind=8),    parameter ::   defepsproba   =  1.00D-14
      real(kind=8),    parameter ::   defepsneg     = -1.00D-6
      real(kind=8),    parameter ::   defepseigen   = 1.00D-4
      real(kind=8),    parameter ::   defsmalleigen = 1.00D-10
      integer(kind=4), parameter ::   maxbondef     =  1
      integer(kind=4), parameter ::   maxmaxb       =  4
      integer(kind=4), parameter ::   mxcentdef     =  6
      integer(kind=4), parameter ::   mxdafhdef     =  200
      integer(kind=4), parameter ::   maxfno        =  100
c
c=======================================================================
c
c     S T A R T
c
c=======================================================================
c
      call timestamp (line)
      call timer (0,ipid,'_edf      ',-1)
      call timer (1,ipid,'_edf      ',-1)
c
c.....factorials
c
      fact(+0)  = 1.0d+00
      facd(-1)  = 1.0d+00
      facd(+0)  = 1.0d+00
      do i=1,maxfac
        fact(i) = dble(i)*fact(i-1)
        facd(i) = dble(i)*facd(i-2)
      enddo
c
c.....Determine the size occupied by a real*8 variable in a direct
c     access file. The routine lengrec () return in the common /sizew/ 
c     the value of the variable RecLength. This variable is then used 
c     in all the routines opening direct access files.
c
      call lengrec ()
c
c.....Logical units. 
c
c     udat   = wfn file 
c     stdout = standard output
c     stdin  = standard input
c     stderr = standard error 
c     lu18   = Atomic Overlap matrix (AOM) file
c     ifilc  = Coefficients of wave function expansion in Slater dets.
c     ifilmo = Idem with Slater dets built in with localized MOs.
c
      udat   =  1
      udatnw =  2
      stdout =  6
      stdin  =  5
      stderr =  0
      lu18   = 18
      ifilc  = 57
      ifilmo = 71
      lw     = stdout
c
c-----------------------------------------------------------------------
c
c.....Logical variables (and their default (D) values)
c
c--Variable--D--Meaning-------------------------------------------------
c
c  skipcan   F  Exact probabilities using the input MOs are computed
c
c  skiplocx  F  Approximate probabilities using Localized MOs are computed
c
c  skiploc   T  Exact probabilities using Localized MOs are not computed
c
c  localize  F  Isopycnic localized MOs of the natural MOs are not computed
c
c  canloc    F  Isopycnic localized Canonical MOs are not computed
c
c  chemiloc  F  Chemical localization is not carried out
c
c  dafhloc   F  Iterative DAFH localization is not carried out
c
c  otherdafh F  The second method of succesive DAFH localizations is not
c               carried out.
c
c  openloc   F  Iterative localization based on open systems is not 
c               carried out
c
c  otheropen F  The second method of iterative localization based on 
c               open systems is neither carried out
c
c  aomloc    F  Iterative AOM localization is not carried out
c
c
c  short     T  Short output requested
c
c  coreloc   F  Core orbitals are not automatically constructed 
c               from Localized MOs
c
c  corecan   F  Core orbitals are not automatically constructed 
c               from canonical MOs
c
c  qtaim     T  QTAIM space partitioning in the computation of the AOM
c
c  mulliken  T  Space partitioning according to the Mulliken population 
c               analysis in the computation of the AOM
c
c  mindef    T  Space partitioning according to the Minimally Deformed 
c               Atoms criterion of Fernandez Rico et al. in the 
c               computation of the AOM()
c
c  mindefrho T  Same as above but using weights 'w_A=rho_A/rho' 
c               that always give S_ii^A>0
c
c  netrho    T  Same as above but using weights 'w_A=rho_A/rho' from 
c               net densities.
c
c  promrho   T  Same as above but using weights 'w_A=rho_A/rho' from 
c               Koga atomic densities 'rho_A' and the promolecular 
c               density 'rho' obtained as 'SUM rho_A'
c
c
c  hesel     T  In a chemical localization the weights defined by 
c               Andreas Heselmann are used. 
c               (A. Heselmann, JCTC 12, 2720-2741 (2016).)
c
c
c  becke     T  Space partitioning according to the  Becke prescription 
c               in the computation of the AOM
c
c  wfncut    F  When epswfn is set to a number greater than 0.0 with 
c               the EPSWFN order, wfncut is set to .true.
c
c
c  lowdin    T  Space partitioning according to the Lowdin population 
c               analysis in the computation of the AOM
c
c  orderp    T  Spinless probabilities are orderer by decreasing value 
c               when its number is smaller than NSTACK. This parameter 
c               is currently set to 4000. 
c
c
c  locore    F  CORE MOs from a multideterminant wavefunction are not 
c               localized before starting the calculations.
c
c
c  memedf    T  For closed-shell systems and multiDet WFNs the ALPHA 
c               and BETA probabilities are saved in arrays instead of
c               being stored in temporary files. This is the default
c               and in this case the routine 'mcalcedf.f' instead of 
c               'calcedf.f' is used. This option can be changed with 
c               the keyword NOMEM.
c
c  rohfloc   F  Localize the ALPHA and BETA spin-orbitals of a ROHF 
c               WFN separately and stop.
c
c  uhfloc    F  Localize the ALPHA and BETA spin-orbitals of a UHF 
c               WFN separately and stop.
c
c
c  opensys   F  Performs an open quantum system (OQS) analysis and stops.
c
c  oqseos    F  Performs an open quantum system (OQS) analysis, computes
c               the effective oxidation state (EOS) of each fragment, and 
c               stops. Actually, it is very similar to the above 'opensys'
c               order, but more recent.
c
c  dolewis   F  Performs a search of Lewis structures based on fragments
c               instead of atoms. When the keyword 'OQSLEWIS' is used in
c               the input, this variable is changed to .true., the routine
c               'lewisoqs.f' is called, and the program stops. The routine
c               lewisoqs.f tries to perform a search of Lewis structures
c               based on the user defined fragments by means of the 
c               'NFRAG' order.
c
c  dofno     F  Performd and Fragment Natural Orbitals (FNO) analysis 
c               and stops
c
c  dodirprob F  Computes by a pure brute force method the probabilities
c               of the sub-set of RSRSs with electron populations between
c               minpopul(i) and maxpopul(i) values. At this time, it only 
c               works for SDWs.
c
c  ditwocent  F The DI between all pairs of atoms is not computed. 
c               With the TWOCENDI order 'ditwocent' is set to .true.
c
c
c  defragoqs  F The values of 'nfrag', 'nwithin(:)', and ifrag(:,:)
c               are taken equal to the values of 'ngroup','nfugrp(:)'
c               and 'ifugrp(:,:)', respectively. To define other values
c               for 'nfrag', 'nwithin(:)', and ifrag(:,:), an ordre
c               'NFRAG' must be used in the input file. When this is done
c               defragoqs is set to .true.
c
c-----------------------------------------------------------------------
c
      largwr    = .false.
      skipcan   = .false.
      skiploc   = .true.
      skiplocx  = .false.
      localize  = .false.
      canloc    = .false.
      chemiloc  = .false.
      dafhloc   = .false.
      otherdafh = .false.
      openloc   = .false.
      otheropen = .false.
      aomloc    = .false.
      short     = .true.
      wfncut    = .false.
      coreloc   = .false.
      corecan   = .false.
      allcloc   = .false.
      orderp    = .true.
      locore    = .false.
      memedf    = .true.
      ditwocent = .false.
      nindep    = .false.
      bonding   = .false.
      maxwfn    = .false.
      aomnorm   = .false.
      recur     = .false.
      sprimok   = .false.
      ngless    = .false.
      opensys   = .false.
      oqseos    = .false.
      dolewis  = .false.
      ok1rdm    = .false.
      dofno     = .false.
      doentropy = .false.
      rohfloc   = .false.
      uhfloc    = .false.
      dodirprob = .false.
      defragoqs = .false.
      salvanu   = .false.
      scuseria  = .false.
      ovcritloc = highover
      ovcritcan = highover
      randommin = ranmindef
      randommax = ranmaxdef 

      pow       = 3     ! Default K value in the iterative Becke scheme
      nrad      = 100   ! Default number of radial points in  "     "
      nang      = 434   ! Default number of Lebedev points    "     "
      irmesh    = 1     ! Default mapping of the radial coordinate r
      rhopow    = 1.0D0 ! Default power of each rho_A when computing w_A
c                         as w_A=rho_A^rhopow/SUM(rho_A^rhopow) 
c                         in the PROMRHO partition.
c
c.....tolsym = Epsilon parameter used in symmetry modules that measures 
c     the tolerance for identity of points,matrices,...
c.....epssym = Parameter that determines if molecular orbitals in transor.f 
c     are considered linearly independent or not.
c.....epsorb = Parameter that determines if molecular orbitals in transor.f
c     break symmetry or not.
c
      tolsym      = 1.0d-6
      epssym      = 1.0d-6
      epsorb      = 1.0d-4
      lpointgroup = .false.
      cpointgroup = 'non'
      nosym       = .false.
      dosymmetry  = .false.
      modelsym    = 1
      alphahesel  = defahesel
      betahesel   = defbhesel
      damps       = defdamps
      covx        = defcovx
      epsneg      = defepsneg
      epseigen    = defepseigen
      smalleigen  = defsmalleigen
      maxbond     = maxbondef
      skiph       = .false.
      dafhselect  = .false.
      mxcent      = mxcentdef
      critov(1:maxcritov)   = critovdef
      ixcent(1:maxcritov)   = .true.
      ndafhsel(1:maxcritov) = 0
      mxdafh           = mxdafhdef
      nfno = 0  ! Number of fragments for which FNOs will be computed.
c
c.....read names of the data files from the standard unit.
c
      write (stdout,1000) line(1:leng(line))
c
c-----Open a file to write the input of this run.
c
      linp = 3
      exfil   = .true.
      fileinp = 'the_input_file'
      do while (exfil)
        inquire (unit=linp,opened=exfil)
        if (.not.exfil) then
          open (unit=linp,file=trim(fileinp))
        else
          linp=linp+1
        endif
      enddo

      do
        read (stdin,'(a)',end=1717) line
        write (linp,'(a)') trim(line)
      enddo
 1717 rewind (stdin)
c
c-----Rewind the file containing the input
c
*     write (stdout,'(a,/,a,/,a,/,a)') ' ',' ',' ',' # '
      write (stdout,'(a)') ' # ++++++++++++++++++++++++++++++++++++++++'
      write (stdout,'(a)') ' # The Input file is'
      write (stdout,'(a)') ' # ++++++++++++++++++++++++++++++++++++++++'
      rewind (linp)
      do
        read (linp,'(a)',end=1718) line
        write (stdout,'(a)') '   '//trim(line)
      enddo
 1718 write (stdout,'(a)') ' # ++++++++++++++++++++++++++++++++++++++++'
      write (stdout,'(a)') ' # End of the Input file'
      write (stdout,'(a)') ' # ++++++++++++++++++++++++++++++++++++++++'
*     write (stdout,'(a)') ' # '
      write (stdout,'(a,/,a,/,a,/,a)') ' # ',' ',' ',' '
      close (linp,status='delete')
c
c-----Start the true reading of the input
c
      read (stdin,'(a)') line
      line=uppcase(line)
      lp=1
      ok=setint(ioverlap,line,lp)
      if (.not.ok) then
        write (stderr,321) 
        stop
      endif
c
      read (stdin,'(a)') filedat
      filedat  = filedat(1:leng(filedat))
      typepart = filedat(1:leng(filedat))
      lp = 1
      ok = setword(word,typepart,lp)
      if (.not.ok) then
        stop ' # edf.f: Second record of input file is empty'
      else
        line = typepart(lp:)
      endif
c
      mulliken  = .false.
      qtaim     = .false.
      mindef    = .false.
      mindefrho = .false.
      netrho    = .false.
      promrho   = .false.
      lowdin    = .false.
      hesel     = .false.
      becke     = .false.
      if (uppcase(typepart(1:8)).eq.'MULLIKEN') then
        mulliken = .true.
      elseif (uppcase(typepart(1:6)).eq.'LOWDIN') then
        lowdin   = .true.
      elseif (uppcase(typepart(1:6)).eq.'MINDEF'.and.
     &        uppcase(typepart(1:9)).ne.'MINDEFRHO') then
        mindef   = .true.
      elseif (uppcase(typepart(1:9)).eq.'MINDEFRHO' .or.
     &        uppcase(typepart(1:5)).eq.'BECKE'     .or.
     &        uppcase(typepart(1:6)).eq.'NETRHO'    .or.
     &        uppcase(typepart(1:7)).eq.'PROMRHO'   .or.
     &        uppcase(typepart(1:9)).eq.'HESELMANN' ) then
        if (uppcase(typepart(1:9)).eq.'MINDEFRHO') mindefrho = .true.
        if (uppcase(typepart(1:5)).eq.'BECKE')     becke     = .true.
        if (uppcase(typepart(1:6)).eq.'NETRHO')    netrho    = .true.
        if (uppcase(typepart(1:7)).eq.'PROMRHO')   promrho   = .true.
        if (uppcase(typepart(1:9)).eq.'HESELMANN') hesel     = .true.
        lp       = 1
        ok       = setint(nang,line,lp)
        if (ok) then
          nang=min(abs(nang),5810)
          nang=max(nang,6)
          ok1 = setint(nrad,line,lp)
          if (ok1) then
            nrad=min(abs(nrad),5000)
            nrad=max(nrad,10)
            ok2 = setint(irmesh,line,lp)
            if (ok2) then
              irmesh=min(abs(irmesh),4)
              irmesh=max(irmesh,1)
              ok3 = setint(pow,line,lp)
              if (ok3) then
               pow=min(abs(pow),9)
               pow=max(pow,2)
               ok4 = setdble(rhopow,line,lp)
               if (ok4) then
                 rhopow = max(0.25D0,rhopow)
                 rhopow = min(2D0,rhopow)
                 ok5 = setdble(alphahesel,line,lp)
                 if (ok5) then
                   alphahesel=abs(alphahesel)
                   ok6 = setdble(betahesel,line,lp)
                   if (ok6) betahesel=abs(betahesel)
                 endif
               endif
              endif
            endif
          endif
        endif
      else
        qtaim    = .true.
      endif
c
      if (mulliken) then
        write (stdout,33) 'MULLIKEN atoms'
      elseif (lowdin) then
        write (stdout,33) 'LOWDIN atoms'
      elseif (mindef) then
        write (stdout,33) 'MINDEF atoms'
      elseif (mindefrho) then
        write (stdout,33) 'MINDEFRHO atoms'
      elseif (netrho) then
        write (stdout,33) 'NETRHO atoms'
      elseif (promrho) then
        write (stdout,33) 'PROMRHO atoms'
      elseif (hesel) then
        write (stdout,33) 'HESELMANN atoms'
      elseif (becke) then
        write (stdout,33) 'BECKE atoms'
      else
        if (ioverlap.ne.-2) then
          write (stdout,33) filedat(1:leng(filedat))
          open (lu18,file=filedat,status='old',iostat=ier)
          if (ier.ne.0) then
            write (stderr,*) '# edf.f: Error opening AOM file'
            stop
          endif
        endif
      endif
c
c-----Here an IF THEN ...ELSE only for non complex AOM.
c
      if (ioverlap.ne.-1.and.ioverlap.ne.-2.and.ioverlap.ne.-11) then
        read (stdin,'(a)') wfnfile
        wfnfile=wfnfile(1:leng(wfnfile))
        open (udat,file=wfnfile,status='old')
        write (stdout,34) wfnfile(1:leng(wfnfile))
c
c.......Read wave function.
c
        call rdwfn (udat,stdout,stderr,ifilc,mal,mbe,wfnfile,cicoef)
c
c-------Probabilities are not possible with CCIQA WFNs
c
        if (cciqa) then
          if ((.not.skipcan).or.(.not.skiploc).or.(skiplocx)) then
            write (stderr,1117) 
            skipcan  = .true.
            skiploc  = .true.
            skiplocx = .true.
          endif
        endif
c
c.......Print coordinates of all the atoms and definition of the groups.
c
        write (stdout,501) 
        write (stdout,503)
        do i=1,ncent
          write (stdout,103) atnam(i)(1:4),i,(xyz(i,k),k=1,3),charge(i)
        enddo
        write (stdout,21) nel,nmo
        if (ndets.gt.1) then
          write (stdout,361) ndets,nelact,ncore,nact
        else
          ncore=0
          write (stdout,36) ndets
          if (rhf)      write (stdout,363) 'RHF '
          if (uhf)      write (stdout,363) 'UHF'
          if (rohf)     write (stdout,363) 'ROHF'
        endif
        if (qtaim.and.largwr) then
          if (ioverlap.eq.0) write (stdout,55) 'PROMOLDEN'
          if (ioverlap.ne.0) write (stdout,55) 'AIMPAC'
          if (ioverlap.eq.2) write (stdout,55) 'TOPMODELF'
          if (ioverlap.eq.3) write (stdout,55) 'TOPMODAIM'
          if (ioverlap.eq.4) write (stdout,55) 'AIMALL'
        endif
        cicoef=cicoef(1:leng(cicoef))
c
c.......Direct access file of wave function coefficients is open.
c
        if ((.not.cciqa).and.(.not.onlynat)) then
          open (ifilc,file=cicoef,access='direct',
     &    recl=RecLength*(nelact+1),form='unformatted')
        endif
c
c.......Read Atomic overlap matrix (AOM). 
c
c       TOPMOD. The number of basins is different to the number of atoms.
c       The file containing the integrals is written in double precision. 
c       Read it and convert it to quadruple.
c
        if (ioverlap.eq.2.or.ioverlap.eq.3) then
          read (lu18,*) nbasins, nmo, natoms
c
c---------This version of EDF only allows 'nbasins' and 'natoms' to both 
c         be equal to 'ncent', the variable read from the wave function 
c         file. In case this is not so, it is necessary to add ghost atoms
c         in this file.
c
          if (nbasins.ne.ncent .or. natoms.ne.ncent) then
            write (stderr,76)
            stop
          endif
        endif
c
        allocate (aom(ncent,nmo,nmo))
        allocate (sgallc(ncent,nmo,nmo))
        allocate (iaom(ncent))
        if (qtaim) then
          if (ioverlap.lt.2.or.ioverlap.eq.4) then
            iaom(1:ncent)=0
            do i=1,ncent
              read (lu18,*) iaom(i)
              if (iaom(i).gt.ncent) write (stderr,17) iaom(i),ncent
              if (iaom(i).gt.ncent) stop
              if (largwr) write (stdout,16) iaom(i)
              ii=iaom(i)
              if (ioverlap.eq.0) then
                read (lu18,80) ((aom(ii,m,j),m=1,j),j=1,nmo) ! PROMOLDEN
              elseif (ioverlap.eq.4) then
                do m=1,nmo
                  read (lu18,81) (aom(iaom(i),j,m),j=1,m)     ! AIMALL
                enddo
              else
                do m=1,nmo
                  read (lu18,82) (aom(iaom(i),j,m),j=1,m)     ! AIMPAC
                enddo
              endif
            enddo
c
c...........topmod
c
          elseif (ioverlap.eq.2.or.ioverlap.eq.3) then    !topmod
            read (lu18,'(a)') line
            lp=1
            do i=1,ncent
              ok = setword(namebas(i),line,lp)
              if (.not.ok) then
                write (stderr,*) '# edf.f: Error reading namebas(',i,')'
                stop
              endif
            enddo
            if (ioverlap.eq.2) write (stdout,1133) 'ELF'
            if (ioverlap.eq.3) write (stdout,1133) 'AIM'
            write (stdout, *) '# Number of elf basins & names', ncent
            do i=1,ncent
              write (stdout,'(a)') ' # '//trim(namebas(i))
            enddo
            allocate (aomtopmod(ncent,nmo,nmo))
            do i=1,ncent
              iaom(i)=i
              read (lu18,80) ((aomtopmod(iaom(i),j,m),m=j,nmo),j=1,nmo)
            enddo
            aom=aomtopmod
            deallocate (aomtopmod)
          endif
          close (lu18)
        else
CC          call timer (2,ifuzan,'_fuzzy    ',-1)
CCc
CCc---------Orthogonalize exactly canonical and natural MOs using OFMO
CCc
CC          write (stdout,10)
CC          allocate (eta(nmo,nmo))
CC          call ofmortho (coef,eta,nmo,nprims,warno,stdout)
CC          if (warno) then
CC             stop ' # edf.f: Singular matrix in ofmortho.f'
CC          else
CC             write (stdout,*) '# OFMO orthogonalization ends Ok!'
CC          endif
CC          deallocate (eta)
CC          forall (i=1:ncent) iaom(i)=i
CC          if (mulliken) then
CC            call aomulliken (aom,0,.false.,stdout,wfnfile)
CC          elseif (mindef) then
CC            call aomindef (aom,0,.false.,stdout,wfnfile)
CC          elseif (mindefrho.or.netrho.or.promrho.or.hesel) then
CC            if (mindefrho) irho=1
CC            if (netrho)    irho=2
CC            if (promrho)   irho=3
CC            if (hesel)     irho=4
CC            allocate (cc(nmo+nmo,nprims))
CC            cc=coef
CC            call aomindefrho (aom,alphahesel,betahesel,
CC     &           irho,nrad,nang,irmesh,rhopow,.false.,cc,stdout,wfnfile)
CC            deallocate (cc)
CC          elseif (lowdin) then
CC            call aomlowdin (aom,0,.false.,stdout,wfnfile)
CC          elseif (becke) then
CC            allocate (cc(nmo,nprims))
CC            cc(1:nmo,1:nprims)=coef(nmo+1:nmo+nmo,1:nprims)
CC            call aombecke 
CC     &        (aom,nrad,nang,pow,irmesh,.true.,cc,stdout,wfnfile)
CC            deallocate (cc)
CC          endif
CC          call timer (4,ifuzan,'_fuzzy    ',-1)
        endif
c
c.......Test if the AOM has been read in or computed for all the basins.
c
        if (qtaim) then
          notaom=.false.
          do i=1,ncent
            j=1
            okaom=.false.
            do while (j.le.ncent.and.(.not.okaom))
              if (iaom(j).eq.i) okaom=.true.
              j=j+1
            enddo
            if (.not.okaom) then
              notaom=.true.
              write (stderr,15) i
              stop
            endif
          enddo
        endif
        if (notaom) stop ' # edf.f: !! Bad definition of the AOM'
c
c.......Fill in the lower part of the aom() array.
c
        if (qtaim) then
          do m=1,nmo
            do j=1,m
              aom(1:ncent,m,j)=aom(1:ncent,j,m)
            enddo
          enddo
        endif
        deallocate (iaom)
      elseif (ioverlap.eq.-1) then
c
c.......COMPLEX AOM CASE
c
        call edfx (stdin,stdout,lu18)
        call timer (5,ipid,'_edf      ',stdout)
        call timestamp (line)
        write (stdout,2000) line(1:leng(line))
        stop 'edf: End of Run Ok !!!'
      elseif (ioverlap.eq.-2) then
c
c.......COMPLEX AOM CASE with AOM integrals from the critic2 program
c
        write (stdout,330) trim(filedat)
        open (unit   = lu18,
     &        file   = trim(filedat),
     &        status = 'old',
     &        form   = 'UNFORMATTED',
     &        iostat = ier)
        if (ier.ne.0) then
          write (stderr,*) '# edf.f: Error opening CRITIC2  AOM file'
        endif
        call edfcrit (stdin,stdout,lu18)
        call timer (5,ipid,'_edf      ',stdout)
        call timestamp (line)
        write (stdout,2000) trim(line)
        stop 'edf: End of Run Ok !!!'
      elseif (ioverlap.eq.-11) then
        call edfstd (stdin,stdout,stderr,lu18)
        call timer (5,ipid,'_edf      ',stdout)
        call timestamp (line)
        write (stdout,2000) line(1:leng(line))
        stop 'edf: End of Run Ok !!!'
      endif
c 
c.....Initizalize some variables
c
      epsdet     = zero
      epsproba   = defepsproba
      tolaom     = deftolaom
      epswfn     = defepswfn
      epswfnloc  = defepswfn
      critics    = defcritics
      probcut    = zero
      nactdet    = ndets
      ngroup     = izero
      mocore     = izero
      moval      = nmo
      nrsrs      = izero
      allocate (ival(nmo))
      allocate (icore(nmo))
      allocate (icorea(nmo))
      allocate (icoreb(nmo))
      forall (i=1:moval) ival(i)=i
c
c-----------------------------------------------------------------------
c
c       BEGIN READING PARTICULAR CONDITIONS OF THE RUN
c
c-----------------------------------------------------------------------
c
      goon=.true.
 20   do while (goon) 
        read (stdin,'(a)',end=110) line
        line=uppcase(line)
        lp=1
        ok=setword(word,line,lp)
        lword=leng(word)
        if (lword.eq.0.or.word(1:1).eq.'#') then
        elseif (word(1:lword).eq.'END') then
          goon=.false.
        elseif (word(1:lword).eq.'PRSRS') then
          nrsrs = nrsrs + 1
          if (nrsrs.le.numrsrs) then
             nres(nrsrs) = line(lp:)
          else
            nrsrs=nrsrs-1
            write (stderr,*) '# edf.f: Too many PRSRS orders'
            write (stderr,*) '# edf.f: Maximum value is ',numrsrs
          endif
          skipcan  = .true.
          skiploc  = .true.
          skiplocx = .true.
        elseif (word(1:lword).eq.'LARGE') then
          short = .false.
          largwr = .true.
        elseif (word(1:lword).eq.'VERBOSE') then
          short = .false.
          largwr = .true.
        elseif (word(1:lword).eq.'BIGOUTPUT') then
          short = .false.
          largwr = .true.
        elseif (word(1:lword).eq.'RECUR') then
          recur = .true.
        elseif (word(1:lword).eq.'AOMNORM') then
          if (qtaim) then
            line=line(lp:)
            call aomnorma 
     &       (line,aom,tolaom,stdout,stderr,lu18,aomnorm,filedat)
          endif
        elseif (word(1:6).eq.'FULLOC' .or. 
     &          word(1:6).eq.'CANLOC' .or.
     &          word(1:6).eq.'AOMLOC' .or.
     &          word(1:7).eq.'DAFHLOC' .or.
     &          word(1:7).eq.'CHEMLOC' .or.
     &          word(1:7).eq.'OPENLOC') then
          if (word(1:6).eq.'FULLOC')  localize = .true.
          if (word(1:6).eq.'CANLOC')  canloc   = .true.
          if (word(1:6).eq.'AOMLOC')  aomloc   = .true.
          if (word(1:7).eq.'DAFHLOC') dafhloc  = .true.
          if (word(1:7).eq.'CHEMLOC') chemiloc = .true.
          if (word(1:7).eq.'OPENLOC') openloc  = .true.
          ok = setword(word,line,lp)
          if (ok) then
            allcloc = .false.
            if (word(1:3).eq.'ALL') allcloc=.true.
          endif
        elseif (word(1:lword).eq.'ROHFLOC') then
          rohfloc = .true.
        elseif (word(1:lword).eq.'UHFLOC') then
          uhfloc = .true.
        elseif (word(1:lword).eq.'DOENTROPY') then
          doentropy = .true.
        elseif (word(1:lword).eq.'DIRPROB') then
          dodirprob = .true.
        elseif (word(1:lword).eq.'DODAFH' .or. 
     &          word(1:lword).eq.'OTHEROPEN') then
           line = line(lp:) 
           call allocate_space_for_dafh (mxdafh,ncent)
           call dafhdrv (stdin,stdout,mxdafh)
           if (word(1:lword).eq.'DODAFH')    otherdafh = .true.
           if (word(1:lword).eq.'OTHEROPEN') otheropen = .true.
           allcloc = .true.
        elseif (word(1:lword).eq.'FNO') then
           if (.not.allocated(inside)) then
             allocate (inside(ncent,maxfno))
             allocate (ninside(maxfno))
             dofno=.true.
           endif
           nfno=nfno+1
           if (nfno.gt.maxfno) then
             stop ' # edf.f: Increase the value of MAXFNO parameter'
           endif
c
c..........Define a fragment to compute its Fragment Natural Orbitals
c
           ninfrag = 0
           do while (setint(iatomo,line,lp))
             if (iatomo.gt.0.and.iatomo.le.ncent) then
               do j=1,ninfrag
                 if (iatomo.eq.inside(j,nfno)) goto 11
               enddo
               ninfrag = ninfrag+1
               inside(ninfrag,nfno)=iatomo
 11            continue
             endif
           enddo
           ninside(nfno)=ninfrag
        elseif (word(1:lword).eq.'OPENSYS') then
          opensys = .true.
        elseif (word(1:lword).eq.'OQSEOS') then
          oqseos = .true.
        elseif (word(1:lword).eq.'OQSLEWIS') then
          dolewis = .true.
        elseif (word(1:lword).eq.'DORDM') then
          ok1rdm = .true.
        elseif (word(1:lword).eq.'NOPROB') then
          skipcan = .true.
        elseif (word(1:lword).eq.'PROBLOC') then
          skiploc = .false.
        elseif (word(1:lword).eq.'NOPROBX') then
          skiplocx = .true.
        elseif (word(1:lword).eq.'PROBCUT') then
          line=line(lp:)
          read (line,*) probcut
        elseif (word(1:lword).eq.'EPSDET') then
          line=line(lp:)
          read (line,*) epsdet
          epsdet=abs(epsdet)
        elseif (word(1:lword).eq.'RHOPOW') then
          line=line(lp:)
          lp=1
          ok=setdble(rhopow,line,lp)
          if (ok) then
            rhopow = abs(rhopow)
            rhopow = max(0.25D0,rhopow)
            rhopow = min(2D0,rhopow)
          endif
        elseif (word(1:lword).eq.'CRITICS') then
          line=line(lp:)
          lp=1
          ok=setdble(critics,line,lp)
          if (ok) then
            critics=abs(critics)
            critics=min(critics,1.0D-1)
            critics=max(critics,1.0D-6)
          endif
        elseif (word(1:lword).eq.'EPSNEG') then
          line = line(lp:)
          lp=1
          ok = setdble(epsneg,line,lp)
          if (ok) epsneg=-abs(epsneg)
        elseif (word(1:lword).eq.'EPSEIGEN') then
          line = line(lp:)
          lp=1
          ok = setdble(epseigen,line,lp)
          if (ok) epseigen=abs(epseigen)
        elseif (word(1:lword).eq.'SMALLEIGEN') then
          line = line(lp:)
          lp=1
          ok = setdble(smalleigen,line,lp)
          if (ok) smalleigen=abs(smalleigen)
        elseif (word(1:lword).eq.'CRITOV') then
          line=line(lp:)
          lp=1
          lpold=lp
          ok=setint(ncritov,line,lp)
          if (ok) then
            if (ncritov.gt.0.and.ncritov.lt.maxcritov) then
              ok=setdble(critval,line,lp)
              if (ok) critov(ncritov)=abs(critval)
            endif
          else
            ok=setdble(critval,line,lpold)
            if (ok) critov(1:maxcritov)=critval
          endif
        elseif (word(1:lword).eq.'CRITOV1') then
          line=line(lp:)
          lp=1
          ok=setdble(critov(1),line,lp)
          if (ok) then
            critov(1)=abs(critov(1))
            if (critov(1).gt.0.995d0) critov(1)=0.995d0
            if (critov(1).lt.0.600d0) critov(1)=0.600d0
          endif
        elseif (word(1:lword).eq.'CRITOV2') then
          line=line(lp:)
          lp=1
          ok=setdble(critov(2),line,lp)
          if (ok) then
            critov(2)=abs(critov(2))
            if (critov(2).gt.0.995d0) critov(2)=0.995d0
            if (critov(2).lt.0.600d0) critov(2)=0.600d0
          endif
        elseif (word(1:lword).eq.'DAMPS') then
          line=line(lp:)
          lp=1
          ok=setdble(damps,line,lp)
          if (ok) then
            damps=abs(damps)
            if (damps.gt.0.999d0) damps=0.999d0
            if (damps.lt.0.800d0) damps=0.800d0
          endif
        elseif (word(1:lword).eq.'MXCENT') then
          line=line(lp:)
          lp=1
          ok=setint(mxcent,line,lp)
          if (ok) then
            line=line(lp:)
            lp=1
            mxcent=abs(mxcent)
            if (mxcent.gt.maxcritov) mxcent=maxcritov 
            if (mxcent.lt.1)  mxcent=1
            ok1=.true.
            do while (ok1) 
              ok1=setint(nn,line,lp)
              if (ok1.and.nn.ge.1.and.nn.le.maxcritov) then
                ixcent(nn) = .false.
              endif
            enddo
          endif
        elseif (word(1:lword).eq.'DAFHSELECT') then
          line=line(lp:)
          lp=1
          if (setint(nn,line,lp)) then
            if (nn.ge.1.and.nn.le.maxcritov) then
              dafhselect(nn) = .true.
              ndafhsel(nn) = ndafhsel(nn) + 1
              inu = ndafhsel(nn)
              nlista=0
              donn: do
                ok=setint(natis,line,lp)
                if (ok) then
                  ok=ok.and.natis.gt.0.and.natis.le.ncent
                  if (ok) then
                    do i=1,nlista
                      if (natis.eq.idafhsel(nn,inu,i)) cycle donn
                    enddo
                    nlista=nlista+1
                    idafhsel(nn,inu,nlista)=natis
                  endif
                else
                  exit donn
                endif
              enddo donn
              if (nlista.ne.nn) then
                stop 'edf.f: Wrong DAFHSELECT order format'
              endif
            endif
          endif
        elseif (word(1:lword).eq.'COVX') then
          line=line(lp:)
          lp=1
          ok=setdble(covx,line,lp)
          if (ok) then
            covx=abs(covx)
            if (covx.gt.5.0d0) covx=5.0d0
            if (covx.lt.0.5d0) covx=0.5d0
          endif
        elseif (word(1:lword).eq.'SKIPH') then
          skiph = .true.
        elseif (word(1:lword).eq.'SALVADORNU') then
          salvanu = .true.
        elseif (word(1:lword).eq.'MAXBOND') then
          line=line(lp:)
          lp=1
          ok=setint(maxbond,line,lp)
          if (ok) then
            maxbond=max(abs(maxbond),maxbondef)
            maxbond=min(maxbond,maxmaxb)
          endif
        elseif (word(1:lword).eq.'ALPHA_HESELMANN') then
          line=line(lp:)
          lp=1
          ok=setdble(alphahesel,line,lp)
          if (ok) alphahesel=abs(alphahesel)
        elseif (word(1:lword).eq.'BETA_HESELMANN') then
          line=line(lp:)
          lp=1
          ok=setdble(betahesel,line,lp)
          if (ok) betahesel=abs(betahesel)
        elseif (word(1:lword).eq.'KPOWER') then
          line=line(lp:)
          lp=1
          ok=setint(pow,line,lp)
          if (ok) then
            pow=min(abs(pow),9)
            pow=max(pow,2)
          endif
        elseif (word(1:lword).eq.'LEBEDEV') then
          line=line(lp:)
          lp=1
          ok=setint(nang,line,lp)
          if (ok) then
            nang=min(abs(nang),5810)
            nang=max(nang,6)
          endif
        elseif (word(1:lword).eq.'NRAD') then
          line=line(lp:)
          lp=1
          ok=setint(nrad,line,lp)
          if (ok) then
            nrad=min(abs(nrad),5000)
            nrad=max(nrad,10)
          endif
        elseif (word(1:lword).eq.'IRMESH') then
          line=line(lp:)
          lp=1
          ok=setint(irmesh,line,lp)
          if (ok) then
            irmesh=min(abs(irmesh),4)
            irmesh=max(irmesh,1)
          endif
        elseif (word(1:lword).eq.'EPSPROBA') then
          line=line(lp:)
          lp=1
          ok=setdble(epsproba,line,lp)
          if (ok) then
            epsproba=min(abs(epsproba),1.0D-1)
          endif
        elseif (word(1:lword).eq.'TOLAOM') then
          line=line(lp:)
          lp=1
          ok=setdble(tolaom,line,lp)
          if (ok) tolaom=max(abs(tolaom),1.0D-3)
        elseif (word(1:lword).eq.'EPSWFN') then
          line=line(lp:)
          read (line,*) epswfnew
          epswfn=max(abs(epswfnew),epswfn)
          if (epswfn.gt.defepswfn) wfncut=.true.
        elseif (word(:lword).eq.'TOLSYM') then
          line=line(lp:)
          read (line,*) tolsym
          tolsym=abs(tolsym)
        elseif (word(:lword).eq.'EPSSYM') then
          line=line(lp:)
          read (line,*) epssym
          epssym=abs(epssym)
        elseif (word(:lword).eq.'EPSORB') then
          line=line(lp:)
          read (line,*) epsorb
          epsorb=abs(epsorb)
        elseif (word(:lword).eq.'DOSYMMETRY') then
          dosymmetry=.true.
          line=line(lp:)
          lp=1
          ok=setint(modelx,line,lp)
          if (ok.and.modelx.ge.1.and.modelx.le.3) then
            modelsym=modelx
          endif
          largwr=.not.short
          call dosym (ncent,nmo,nprims,neq,stdout,modelsym,largwr)
        elseif (word(:lword).eq.'POINTGROUP') then
          line=line(lp:)
          read (line,'(a)') cpointgroup
          cpointgroup=lower(trim(cpointgroup))
          lpointgroup=.true.   
        elseif (word(:lword).eq.'NOSYMMETRY') then
          nosym=.true.
        elseif (word(1:lword).eq.'EPSWFNLOC') then
          line=line(lp:)
          read (line,*) epswfnew
          epswfnloc=max(abs(epswfnew),epswfnloc)
        elseif (word(1:lword).eq.'BONDING') then
          bonding=.true.
        elseif (word(1:lword).eq.'TWOCENDI') then
          ditwocent = .true.
          epsbond=one
          line=line(lp:)
          lp=1
          ok=setdble(epsbond,line,lp)
          if (ok) epsbond=abs(epsbond)
        elseif (word(1:lword).eq.'NDETS') then
          line=line(lp:)
          read (line,*) nactdet
          nactdet=min(iabs(nactdet),ndets)
          if (icorr) write (stdout,28) nactdet
c
c---------Orthogonalize exactly canonical and natural MOs using OFMO
c
        elseif (word(1:lword).eq.'OFMO') then
          write (stdout,10)
          allocate (eta(nmo,nmo))
          call ofmortho (coef,eta,nmo,nprims,warno,stdout)
          if (warno) stop '# edf.f: Singular matrix in ofmortho.f'
          write (stdout,*) '# OFMO orthogonalization ends Ok!'
          deallocate (eta)
        elseif (word(1:lword).eq.'LOCORE') then
           if (ncore.gt.0) locore=.true.
        elseif (word(1:lword).eq.'LOCSOME') then
           if (qtaim) then
             line=line(lp:)
             call locsome (aom,line,stdout,.not.short)
           endif
        elseif (word(1:lword).eq.'NOMEM') then
           memedf=.false.
        elseif (word(1:lword).eq.'CGEDF') then
           if (qtaim) then
             line=line(lp:)
             if (.not.allocated(nfugrp)) 
     &         stop '# edf.f: CGEDF order: nfugrp() is not allocated'
             if (.not.allocated(ifugrp))
     &         stop '# edf.f: CGEDF order: ifugrp() is not allocated'
             call cgedfp (line,wfnfile,
     &        probcut,ngroup,stdout,aom,nfugrp,ifugrp,orderp)
           endif
        elseif (word(1:lword).eq.'RHOCOND') then
           if (qtaim) then
             line=line(lp:)
             if (.not.allocated(nfugrp)) 
     &         stop '# edf.f: RHOCOND order: nfugrp() is not allocated'
             if (.not.allocated(ifugrp))
     &         stop '# edf.f: RHOCOND order: ifugrp() is not allocated'
             call rhocond
     &       (line,wfnfile,udat,ngroup,stdout,aom,nfugrp,ifugrp)
           endif
        elseif (word(1:lword).eq.'ALLRHOC') then
           if (qtaim) then
             if (.not.allocated(nfugrp)) 
     &         stop '# edf.f: ALLRHOC order: nfugrp() is not allocated'
             if (.not.allocated(ifugrp))
     &         stop '# edf.f: ALLRHOC order: ifugrp() is not allocated'
             call allrhoc (wfnfile,udat,ngroup,stdout,aom,nfugrp,ifugrp)
           endif
        elseif (word(1:lword).eq.'RHOALLS') then
           if (qtaim) then
             line=line(lp:)
             if (.not.allocated(nfugrp)) 
     &         stop '# edf.f: RHOALLS order: nfugrp() is not allocated'
             if (.not.allocated(ifugrp))
     &         stop '# edf.f: RHOALLS order: ifugrp() is not allocated'
             call rhoalls (line,wfnfile,
     &        udat,ngroup,stdout,aom,nfugrp,ifugrp)
           endif
        elseif (word(1:lword).eq.'COREMO') then
           coreloc=.true.
           corecan=.true.
           line=line(lp:)
           lp=1
           ok=setdble(ovcritcan,line,lp)
           if (ok) then
             ovcritcan=max(abs(ovcritcan),smallover)
             ovcritloc=max(abs(ovcritcan),smallover)
           endif
        elseif (word(1:lword).eq.'NOORDER') then
          orderp = .false.
        elseif (word(1:lword).eq.'SUPPRESS') then
           if (qtaim) then
             if (ndets.ne.1) then
               line=line(lp:)
               call suppcas (aom,udat,ifilc,epsdet,filedat,
     &              wfnfile,line,stdout)
               call timer (5,ipid,'_edf      ',stdout)
               call timestamp (line)
               write (stdout,2000) trim(line)
               stop 'edf: End of Run Ok !!!'
*              nline='AOM '//trim(line)
*              call suppcas (aom,udat,ifilc,epsdet,filedat,
*    &              wfnfile,nline,stdout)
*              nline='WFN '//trim(line)
*              call suppcas (aom,udat,ifilc,epsdet,filedat,
*    &              wfnfile,nline,stdout)
*              call timer (5,ipid,'_edf      ',stdout)
*              call timestamp (line)
*              write (stdout,2000) trim(line)
*              stop 'edf: End of Run Ok !!!'
             else
              line=line(lp:)
              nline='AOM '//trim(line)
              call suppress (aom,udat,filedat,wfnfile,nline,stdout)
              nline='WFN '//trim(line)
              call suppress (aom,udat,filedat,wfnfile,nline,stdout)
              call timer (5,ipid,'_edf      ',stdout)
              call timestamp (line)
              write (stdout,2000) trim(line)
              stop 'edf: End of Run Ok !!!'
             endif
           endif
        elseif (word(1:lword).eq.'MXDAFH') then
          ok=setint(mxdafh,line,lp)
          if (ok) mxdafh=abs(mxdafh)
        elseif (word(1:lword).eq.'DELETE') then
           if (qtaim) then
             line=line(lp:)
             allocate (aomx(ncent,nmo,nmo))
             aomx(1:ncent,1:nmo,1:nmo)=aom(1:ncent,1:nmo,1:nmo)
             call delmo (aomx,mal,mbe,line,stdout)
             aom(1:ncent,1:nmo,1:nmo)=aomx(1:ncent,1:nmo,1:nmo)
             deallocate (aomx)
           endif
        elseif (word(1:lword).eq.'NFRAG') then
          defragoqs = .true.
          ok=setint(nfrag,line,lp)
          if (ok) then
            nfrag=min(ncent,abs(nfrag))
            allocate (nwithin(nfrag))
            allocate (ifrag(ncent,nfrag))
            do i=1,nfrag
              read (stdin,'(a)') line
              lp=1
              ok=setint(nwithin(i),line,lp)
              if (.not.ok) then
                write (stderr,*) '# edf.f: Error defining fragment',i
                stop
              endif
              do j=1,nwithin(i)
                ok = setint(ifrag(j,i),line,lp)
                if (.not.ok) then
                  write (stderr,*) '# edf.f: Error defining fragment',i
                  stop
                endif
              enddo
c
c.............Test that the group is not empty of atoms.
c
              if (nwithin(i).eq.0) then
                write (stderr,*) '# edf.f: Group ',i,' empty of atoms'
                stop ' # edf.f: GROUP empty of atoms.'
              endif
c
c.............Test that each atom only belongs to a group.
c
              if (i.gt.1) then
                 do k=1,i-1
                   do m=1,nwithin(k)
                     do l=1,nwithin(i)
                       if (ifrag(l,i).eq.ifrag(m,k)) then
                          write (stderr,432) l,i,m,k
                          stop
                       endif
                     enddo
                   enddo
                 enddo
              endif
            enddo
          else
            nfrag=ncent
            allocate (nwithin(nfrag))
            allocate (ifrag(ncent,nfrag))
            forall (i=1:nfrag) nwithin(i)=1
            forall (i=1:nfrag) ifrag(1,i) = i
          endif
        elseif (word(1:lword).eq.'NGROUP' .or.
     &          word(1:lword).eq.'NINDEP') then
          if (word(1:lword).eq.'NINDEP') nindep   = .true.
          ok=setint(ngroup,line,lp)
          if (ok.and.ngroup.lt.0) then
            ngroup=-ngroup
            ngless=.true.
          endif
          ok=ok.and.ngroup.gt.0 
          if (ok) then
            ngroup=abs(ngroup)
c
c...........Maximum number of fragments = number of atoms.
c
            if (ngroup.gt.ncent) write (stderr,383)
            if (ngroup.gt.ncent) stop
            allocate (ifugrp(ncent,ngroup))
            allocate (nfugrp(ngroup))
            allocate (eleca(2,ngroup))
            allocate (elecb(2,ngroup))
            nfugrp  = izero
            ngroupto = ngroup
            if (ngless) ngroupto=ngroup-1
            ntot=0
            do i=1,ngroupto
              read (stdin,'(a)') line
              lp=1
              ok=setint(nfugrp(i),line,lp)
              if (.not.ok) then
                write (stderr,*) '# edf.f: Error defining fragment',i
                stop
              endif
              ntot=ntot+nfugrp(i)
              do j=1,nfugrp(i)
                ok=setint(ifugrp(j,i),line,lp)
                if (.not.ok) then
                  write (stderr,*) '# edf.f: Error defining fragment',i
                  stop
                endif
              enddo
              eleca(1,i)=izero
              elecb(1,i)=izero
              eleca(2,i)=nalpha
              elecb(2,i)=nbeta
c
c.............Test that the group is not empty of atoms.
c
              if (nfugrp(i).eq.0) then
                write (stderr,*) '# edf.f: Group ',i,' empty of atoms'
                stop ' # edf.f: GROUP empty of atoms.'
              endif
c
c.............Test that each atom only belongs to a group.
c
              if (i.gt.1) then
                 do k=1,i-1
                   do m=1,nfugrp(k)
                     do l=1,nfugrp(i)
                       if (ifugrp(l,i).eq.ifugrp(m,k)) then
                         write (stderr,432) l,i,m,k
                         stop
                       endif
                     enddo
                   enddo
                 enddo
              endif
            enddo
            if (ngless) then
               nfugrp(ngroup) = izero
              eleca(1,ngroup) = izero
              elecb(1,ngroup) = izero
              eleca(2,ngroup) = nalpha
              elecb(2,ngroup) = nbeta
              cyclem: do m=1,ncent
                do i=1,ngroup-1
                  do l=1,nfugrp(i)
                    if (ifugrp(l,i).eq.m) cycle cyclem
                  enddo
                enddo
                nfugrp(ngroup)=nfugrp(ngroup) + 1
                ifugrp(nfugrp(ngroup),ngroup) = m
              enddo cyclem
            endif
          endif
c
c         Default for maximum and minimum values total electron populations
c
          allocate (maxpopul(ngroup))
          allocate (minpopul(ngroup))
          maxpopul(:)=nel
          minpopul(:)=0
c
c         Default for maximum and minimum values ALPHA and BETA electron populations
c
          allocate (maxalpha(ngroup))
          allocate (maxbeta(ngroup))
          allocate (minalpha(ngroup))
          allocate (minbeta(ngroup))
          maxalpha(:)= mal
          maxbeta(:) = mbe
          minalpha(:)= 0
          minbeta(:) = 0
        elseif (word(1:lword).eq.'MAXALPHA') then
          if (ngroup.gt.0) then
            line=line(lp:)
            lp=1
            do i=1,ngroup
              ok=setint(maxalpha(i),line,lp)
              if (.not.ok) stop ' # edf.f: Wrong MAXALPHA order'
              maxalpha(i)=min(mal,maxalpha(i))
            enddo
          endif
        elseif (word(1:lword).eq.'MAXBETA') then
          if (ngroup.gt.0) then
            line=line(lp:)
            lp=1
            do i=1,ngroup
              ok=setint(maxbeta(i),line,lp)
              if (.not.ok) stop ' # edf.f: Wrong MAXBETA order'
              maxbeta(i)=min(mbe,maxbeta(i))
            enddo
          endif
        elseif (word(1:lword).eq.'MINALPHA') then
          if (ngroup.gt.0) then
            line=line(lp:)
            lp=1
            do i=1,ngroup
              ok=setint(minalpha(i),line,lp)
              if (.not.ok) stop ' # edf.f: Wrong MINALPHA order'
              minalpha(i)=max(0,minalpha(i))
            enddo
          endif
        elseif (word(1:lword).eq.'MINBETA') then
          if (ngroup.gt.0) then
            line=line(lp:)
            lp=1
            do i=1,ngroup
              ok=setint(minbeta(i),line,lp)
              if (.not.ok) stop ' # edf.f: Wrong MINBETA order'
              minbeta(i)=max(0,minbeta(i))
            enddo
          endif
        elseif (word(1:lword).eq.'MAXPOP') then
            if (ngroup.gt.0) then
              if (.not.allocated(maxpopul)) allocate (maxpopul(ngroup))
              ll=lword
              line=line(lp:)
              read (line,*,iostat=ier) (maxpopul(i),i=1,ngroup)
            else
              write (stderr,*) '# NGROUP=0. Use a previous NGROUP order'
              stop
            endif
*         if (ngroup.gt.0.and.nrsrs.gt.0) then
*           if (.not.allocated(maxpopul)) allocate (maxpopul(ngroup))
*           ll=lword
*           line=line(lp:)
*           read (line,*,iostat=ier) (maxpopul(i),i=1,ngroup)
*         else
*           if (ngroup.le.0) then
*             write (stderr,*) '# NGROUP=0. Use a previous NGROUP order'
*             stop
*           endif
*           if (.not.dodirprob) then
*             if (nrsrs.le.0) then
*               write (stderr,*) '# A previous PRSRS order is needed'
*               stop
*             endif
*           endif
*         endif
        elseif (word(1:lword).eq.'MINPOP') then
            if (ngroup.gt.0) then
              if (.not.allocated(minpopul)) allocate (minpopul(ngroup))
              ll=lword
              line=line(lp:)
              read (line,*,iostat=ier) (minpopul(i),i=1,ngroup)
            else
              write (stderr,*) '# NGROUP=0. Use a previous NGROUP order'
              stop
            endif
*         if (ngroup.gt.0.and.nrsrs.gt.0) then
*           if (.not.allocated(minpopul)) allocate (minpopul(ngroup))
*           ll=lword
*           line=line(lp:)
*           read (line,*,iostat=ier) (minpopul(i),i=1,ngroup)
*         else
*           if (ngroup.le.0) then
*             write (stderr,*) '# NGROUP=0. Use a previous NGROUP order'
*             stop
*           endif
*           if (.not.dodirprob) then
*             if (nrsrs.le.0) then
*               write (stderr,*) '# A previous PRSRS order is needed'
*               stop
*             endif
*           endif
*         endif
        elseif (word(1:lword).eq.'RANDOM') then
          line=line(lp:)
          lp=1
          ok = setdble(ranmin,line,lp)
          if (ok) ok=setdble (ranmax,line,lp)
          randommax=max(ranmin,ranmax)
          randommin=min(ranmin,ranmax)
        elseif (word(1:lword).eq.'MAXWFN') then
          maxwfn=.true.
c
c-----------------------------------------------------------------------
c
c......The following keywords have no sense here but they have in the 
c      'optedf.f' module. Here, DO NOTHING.
c
        elseif (word(1:lword).eq.'ENDBONDING') then
        elseif (word(1:lword).eq.'WEDF')       then
        elseif (word(1:lword).eq.'WPOP')       then
        elseif (word(1:lword).eq.'WDIS')       then
        elseif (word(1:lword).eq.'EPSBOND')    then
        elseif (word(1:lword).eq.'INICVAR')    then
        elseif (word(1:lword).eq.'PRIN')       then
        elseif (word(1:lword).eq.'TYPE')       then
        elseif (word(1:lword).eq.'TYPE3C')     then
        elseif (word(1:lword).eq.'PAIR')       then
        elseif (word(1:lword).eq.'TRIO')       then
        elseif (word(1:lword).eq.'PLUS')       then
        elseif (word(1:lword).eq.'MINUS')      then
c
c.......The following keywords have no sense here but they have in the 
c       'wfnmax.f' module. Here, DO NOTHING.
c
        elseif (word(1:lword).eq.'GTOL')       then
        elseif (word(1:lword).eq.'MAXIT')       then
        elseif (word(1:lword).eq.'MUP')       then
        elseif (word(1:lword).eq.'IPRINT')       then
        elseif (word(1:lword).eq.'CONFIG')       then
c
c-----------------------------------------------------------------------
c
        else
          write (stderr,'(a)') ' # '//word(1:lword)
          write (stderr,*) '# Key word in input file not understood'
        endif
      enddo
 110  continue



c
c-----AOM integrals in fuzzy partitions are computed.
c
      if (.not.qtaim) then
        call timer (2,ifuzan,'_fuzzy    ',-1)
c
c-------Orthogonalize exactly canonical and natural MOs using OFMO
c
        write (stdout,10)
        allocate (eta(nmo,nmo))
        call ofmortho (coef,eta,nmo,nprims,warno,stdout)
        if (warno) then
           stop ' # edf.f: Singular matrix in ofmortho.f'
        else
           write (stdout,*) '# OFMO orthogonalization ends Ok!'
        endif
        deallocate (eta)
        if (mulliken) then
          call aomulliken (aom,0,.false.,stdout,wfnfile)
        elseif (mindef) then
          call aomindef (aom,0,.false.,stdout,wfnfile)
        elseif (mindefrho.or.netrho.or.promrho.or.hesel) then
          if (mindefrho) irho=1
          if (netrho)    irho=2
          if (promrho)   irho=3
          if (hesel)     irho=4
          allocate (cc(nmo+nmo,nprims))
          cc=coef
          call aomindefrho (aom,alphahesel,betahesel,
     &         irho,nrad,nang,irmesh,rhopow,.false.,cc,stdout,wfnfile)
          deallocate (cc)
        elseif (lowdin) then
          call aomlowdin (aom,0,.false.,stdout,wfnfile)
        elseif (becke) then
          allocate (cc(nmo,nprims))
          cc(1:nmo,1:nprims)=coef(nmo+1:nmo+nmo,1:nprims)
          call aombecke (aom,nrad,nang,pow,irmesh,.false.,
     &      cc,stdout,wfnfile)
          deallocate (cc)
        endif
c
c-------Some orders involving the AOM elements need to be perfomed
c       again in case of fuzzy partitions of R^3
c
        rewind (stdin)
        goon=.true.
 201    do while (goon) 
          read (stdin,'(a)',end=111) line
          line = uppcase(line)
          lp=1
          ok = setword(word,line,lp)
          lword=leng(word)
          if (lword.eq.0.or.word(1:1).eq.'#') then
          elseif (word(1:lword).eq.'AOMNORM') then
            line=line(lp:)
            call aomnorma 
     &       (line,aom,tolaom,stdout,stderr,lu18,aomnorm,filedat)
          elseif (word(1:lword).eq.'LOCSOME') then
            line=line(lp:)
            call locsome (aom,line,stdout,.not.short)
          elseif (word(1:lword).eq.'CGEDF') then
            line=line(lp:)
            if (.not.allocated(nfugrp)) 
     &        stop '# edf.f: CGEDF order: nfugrp() is not allocated'
            if (.not.allocated(ifugrp))
     &        stop '# edf.f: CGEDF order: ifugrp() is not allocated'
            call cgedfp (line,wfnfile,
     &      probcut,ngroup,stdout,aom,nfugrp,ifugrp,orderp)
          elseif (word(1:lword).eq.'RHOCOND') then
            line=line(lp:)
            if (.not.allocated(nfugrp)) 
     &        stop '# edf.f: RHOCOND order: nfugrp() is not allocated'
            if (.not.allocated(ifugrp))
     &        stop '# edf.f: RHOCOND order: ifugrp() is not allocated'
            call rhocond
     &      (line,wfnfile,udat,ngroup,stdout,aom,nfugrp,ifugrp)
          elseif (word(1:lword).eq.'ALLRHOC') then
            if (.not.allocated(nfugrp)) 
     &        stop '# edf.f: ALLRHOC order: nfugrp() is not allocated'
            if (.not.allocated(ifugrp))
     &        stop '# edf.f: ALLRHOC order: ifugrp() is not allocated'
            call allrhoc (wfnfile,udat,ngroup,stdout,aom,nfugrp,ifugrp)
          elseif (word(1:lword).eq.'RHOALLS') then
            line=line(lp:)
            if (.not.allocated(nfugrp)) 
     &        stop '# edf.f: RHOALLS order: nfugrp() is not allocated'
            if (.not.allocated(ifugrp))
     &        stop '# edf.f: RHOALLS order: ifugrp() is not allocated'
            call rhoalls (line,wfnfile,
     &       udat,ngroup,stdout,aom,nfugrp,ifugrp)
          elseif (word(1:lword).eq.'SUPPRESS') then
            if (ndets.ne.1) then
              line=line(lp:)
              call suppcas (aom,udat,ifilc,epsdet,filedat,
     &             wfnfile,line,stdout)
              call timer (5,ipid,'_edf      ',stdout)
              call timestamp (line)
              write (stdout,2000) trim(line)
              stop 'edf: End of Run Ok !!!'
            else
             line=line(lp:)
             nline='AOM '//trim(line)
             call suppress (aom,udat,filedat,wfnfile,nline,stdout)
             nline='WFN '//trim(line)
             call suppress (aom,udat,filedat,wfnfile,nline,stdout)
             call timer (5,ipid,'_edf      ',stdout)
             call timestamp (line)
             write (stdout,2000) trim(line)
             stop 'edf: End of Run Ok !!!'
            endif
          elseif (word(1:lword).eq.'DELETE') then
            line=line(lp:)
            allocate (aomx(ncent,nmo,nmo))
            aomx(1:ncent,1:nmo,1:nmo)=aom(1:ncent,1:nmo,1:nmo)
            call delmo (aomx,mal,mbe,line,stdout)
            aom(1:ncent,1:nmo,1:nmo)=aomx(1:ncent,1:nmo,1:nmo)
            deallocate (aomx)
          elseif (word(1:lword).eq.'END') then
              goon=.false.
          else
          endif
        enddo
 111    continue
        call timer (4,ifuzan,'_fuzzy    ',-1)
      endif
c
c-----------------------------------------------------------------------
c       END READING PARTICULAR CONDITIONS OF THE RUN
c-----------------------------------------------------------------------
c
c.....Fill in the AOM() for alpha and beta spin-orbitals, necessary
c     in case of UHF functions.
c 
      call aomspin (ncent,nmo,aom,nalpha,nbeta)
c
c-----Write AOM if large output is active
c
      if (largwr) call writeaom (aom,ncent,nmo,stdout,stderr,uhf)
c
c.....Test again consistency of AOM elements in case AOMNORM keyword
c     has not been read in and the WFN is not an UHF one.
c
      if (.not.uhf) call testaom (aom,tolaom,ncent,nmo,stderr)
c
c.....If none orbital localization is carried out, skip computation
c     of the exact and approximate EDF in terms of localized MOs. 
c
      noloc = .not.localize
      noloc = noloc .and. (.not.canloc    )
      noloc = noloc .and. (.not.chemiloc  )
      noloc = noloc .and. (.not.dafhloc   )
      noloc = noloc .and. (.not.aomloc    )
      noloc = noloc .and. (.not.otherdafh )
      if (noloc) then
        skiploc  = .true.
        skiplocx = .true.
      endif
c
c.....Localize CORE MOs replacing the original MOs of the wavefuncion
c
      if (locore.and.ncore.gt.0) call isopycore (aom,stdout,largwr)
c
c.....No matter whether coreloc is .true. or .false. the
c     CORE aproximations CAN NOT BE USED with multideterminantal
c     wavefunctions made of localized MOs.
c
      if (ndets.ne.1) coreloc = .false.

      if ((.not.skiploc).and.(.not.skiplocx))
     & write (stdout,37) epsdet,epswfn,epswfnloc,probcut,epsproba
c
c.....If ngroup <= 0 ==> each atom is a group with a minimum (maximum)
c     electronic population equal to 0 (nel).
c
      if (ngroup.le.izero) then
        ngroup=ncent 
        if (.not.allocated(ifugrp)) then
          allocate (ifugrp(ncent,ngroup),stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot allocate ifugrp()'
        endif
        if (.not.allocated(nfugrp)) then
          allocate (nfugrp(ngroup),stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot allocate nfugrp()'
        endif
        if (.not.allocated(eleca)) then
          allocate (eleca(2,ngroup),stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot allocate eleca()'
        endif
        if (.not.allocated(elecb)) then
          allocate (elecb(2,ngroup),stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot allocate elecb()'
        endif
        do i=1,ngroup
          nfugrp(i)=1
          ifugrp(1,i)=i
          eleca(1,i)=izero
          elecb(1,i)=izero
          eleca(2,i)=nalpha
          elecb(2,i)=nbeta
        enddo
      endif
c
c.....allocate and compute arrays for group overlap integrals 
c     between canonical MOs.
c
      call allocatesg (ngroup,nmo,nalpha,nbeta)
c
c.....Obtain group overlap integrals.
c
      call computesg (ncent,ngroup,nmo,nalpha,nbeta,aom,nfugrp,ifugrp)
c
c.....If the wavefunction has a single determinant, the canonical
c     and natural MOs are equal and the group overlap integrals
c     are the same. 
c
      if (.not.icorr) sgnat  = sg
      if (     icorr) numdet = nactdet
      if (.not.icorr) numdet = 1
      largwr=.not.short
c
c.....Cut the wavefunctions elliminating coefficients with an absolute 
c     value smaller than epswfn. The NUMDET parameter may possibly
c     differ before and after calling the cutwfn routine.
c
      if (icorr.and.wfncut) then
        numdetprev=numdet
        call cutwfn (epswfn,ifilc,numdet,nelact)
        write (stdout,87) epswfn,numdetprev,numdet
      endif
c
      if (.not.allocated(v1mata)) then
        allocate (v1mata (nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate v1mata()'
      endif
      if (.not.allocated(d1mata)) then
        allocate (d1mata (nmo),    stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate d1mata()'
      endif
c
c-----Compute 1-RDM, update natural MOs and compute AOM between
c     these updated natural MOs.
c
      if (.not.onlynat) then
        if (.not.cciqa) then
          call srhon (epsdet,ifilc,numdet,ncore,nelact,nmo)
        endif
      else
        c1et = 0.0d+00
        forall (i=1:nmo) c1et(i,i)=occ(i)
        ok1rdm = .true.
      endif
      if (largwr.or.ok1rdm) then
        write (lw,*) '#'
        write (lw,*) '# TOTAL 1st-order density matrix'
        write (lw,*) '# (( rho (k,l), k=1,nmo), l=1,k)'
        do i=1,nmo
          write (lw,4444) (c1et(i,j),j=1,i)
        enddo
      endif
c
c.....diagonalize total 1st-order matrix
c
c-----Modification to avoid that c1et is destroyed in the diagonalization
c     by using the auxiliary array eta() to store c1et()
c
      allocate (eta(nmo,nmo))
      eta(1:nmo,1:nmo) = c1et(1:nmo,1:nmo)
c
c-----The following line ....
      if (allocated(work)) deallocate(work)
      allocate (work(max(1,3*nmo-1)),stat=ier)
      if (ier.ne.0) stop 'edf.f: Cannot allocate work() for dsyev'
      call dsyev ('V','U',nmo,eta,nmo,d1mata,work,max(1,3*nmo-1),info)
      if (info.ne.0) write (lw,*) '# edf.f: dsyev failed, info=',info
      deallocate (work)
      v1mata(1:nmo,1:nmo) = eta(1:nmo,1:nmo)
c
c-----replaces the previous jacobi call
*     call jacobi (c1et,nmo,nmo,d1mata,v1mata,nrot)
      deallocate (eta)
c
      if (largwr.or.ok1rdm) then
        write (lw,*) '#'
        write (lw,*) '# First-order density matrix Eigenvalues'
*       write (lw,444) (d1mata(i),i=1,nmo)
        write (lw,4444) (d1mata(i),i=1,nmo)
      endif
c
c.....Discard previously read in natural occupation numbers
c
      forall (i=1:nmo) occ(i) = d1mata(i)
      if (allocated(d1mata)) then
        deallocate(d1mata,stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot deallocate d1mata()'
      endif
      if (largwr) write (lw,*) '# Eigenvectors'
c            
c.....rotate into primitive basis, discarding previous naturals
c       
      do i=1,nmo
        if (largwr) write (lw,'(a,I3)') ' # Eigenvector ',i
        if (largwr) write (lw,4444) (v1mata(j,i),j=1,nmo)
        do j=1,nprims
          tmp = 0d0
          do k=1,nmo
            tmp = tmp + v1mata(k,i) * coef(k+nmo,j)
          enddo
          coef(i,j) = tmp
        enddo
      enddo
      write (lw,*) '# Natural MOs in the WFN file have been replaced'
c
c.....Obtain group overlap integrals between natural MOs
c
      call computesgnat
     &  (ncent,ngroup,nmo,nalpha,nbeta,nfugrp,ifugrp,v1mata) 
c
c-----Localize ALPHA and BETA MOs separately in case of ROFH or UHF WFNs
c
      if (rohfloc.or.uhfloc) then
        if (rohfloc.and.(.not.rohf)) then
          stop 'edf: ROHFLOC order requested but the WFN is not ROHF'
        endif
        if (uhfloc.and.(.not.uhf)) then
          stop 'edf: UHFLOC order requested but the WFN is not ROHF'
        endif
        call locrohf (aomalpha,aombeta,coef,udat,ialpha,ibeta,
     &     ovcritloc,ncent,nmo,nprims,nalpha,nbeta,stdout,
     &     filedat,wfnfile,.false.,largwr)
        call timer (5,ipid,'_edf      ',stdout)
        call timestamp (line)
        write (stdout,2000) line(1:leng(line))
        stop 'edf: End of Run Ok !!!'
      endif
c
c-----Pure brute force calculation of probabilities for a subset of RSRSs
c
      if (.not.onlynat) then
        if (dodirprob) then
          if (ndets.ne.1) then
            stop '# edf.f: DIRPROB order works only for SDWs'
          else
            call dirprob (probcut,nalpha,nbeta,ngroup,stdout,
     &           sgalpha,sgbeta,minpopul,maxpopul)
          endif
          call timer (5,ipid,'_edf      ',stdout)
          call timestamp (line)
          write (stdout,2000) line(1:leng(line))
          stop 'edf: End of Run Ok !!!'
        endif
c
c-------Open systems analysis
c
        if (opensys) then
           woccup = .false.
           call sopenrho (aom,nfugrp,ifugrp,nprims,ncent,nmo,ngroup,
     &       mal,mbe,epsneg,lw,stderr,ifilc,udat,wfnfile,largwr)
           call timer (5,ipid,'_edf      ',stdout)
           call timestamp (line)
           write (stdout,2000) line(1:leng(line))
           stop 'edf: End of Run Ok !!!'
        endif
      endif
c
c-----Open systems analysis & Computation of Effective Oxidation states
c
      if (oqseos.or.dolewis) then

        if (.not.defragoqs) then
          nfrag = ngroup
          allocate (nwithin(nfrag))
          allocate (ifrag(ncent,nfrag))
          nwithin(1:nfrag) = nfugrp(1:ngroup)
          ifrag(1:ncent,1:nfrag) = ifugrp(1:ncent,1:ngroup)
        endif
        if (onlynat) then
          nalpha = nmo
          nbeta = nmo
          c1ea = c1et * 0.5d+00
          c1eb = c1et * 0.5d+00
          forall (i=1:nalpha) ialpha(i) = i
          forall (i=1:nbeta )  ibeta(i) = i
          trzrdma = 0.0d0
          trzrdmb = 0.0d0
          do i=1,nmo
            trzrdma = trzrdma + c1ea(i,i)
            trzrdmb = trzrdmb + c1eb(i,i)
          enddo
          mal = int(anint(trzrdma))
          mbe = int(anint(trzrdmb))
c
c---------If only natural MOs have been read in, the AOM integrals that
c         must be passed to the 'nfnoeos' or 'ulewisoqs' routines are 
c         those between the natural MOs. Otherwise, the AOM integrals
c         must be those between canonical MOs.
c
c---------Compute AOM integrals between natural orbitals
c
          allocate (v1matat(nmo,nmo))
          allocate (aomnat(ncent,nmo,nmo))
          v1matat=transpose(v1mata)
          do i=1,ncent
            aomnat(i,:,:) = matmul(matmul(v1matat,aom(i,:,:)),v1mata)
          enddo
c
c---------aomnat now contains the AOM between natural MOs !!
c
          deallocate (v1matat)
        else
          allocate (aomnat(ncent,nmo,nmo))
          aomnat = aom ! aomnat now contains the AOM between CANMOs !!
        endif
        if (oqseos) then
c
c---------Critical value of the overlap between two normalized in R^3 
c         FNOs that determines whether one considers that they almost
c         the same or not.
c
          critoverlap = 0.95d+00
          if (rhf.or.rohf.or.uhf) then
            call sdweos (aomnat,nwithin,ifrag,mal,mbe,nprims,nmo,
     &      ncent,nfrag,covx,critoverlap,udat,wfnfile,stdout,
     &      largwr)
          else
            call nfnoeos (aomnat,nwithin,ifrag,cciqa,mal,mbe,nprims,nmo,
     &      ncent,nfrag,covx,critoverlap,udat,wfnfile,stdout,stderr,
     &      largwr)
          endif
          call timer (5,ipid,'_edf      ',stdout)
        else
          if (rhf.or.rohf.or.uhf) then
            call sdwlewis (aomnat,mal,mbe,nfrag,nwithin,ifrag,
     &      mxcent,covx,epsneg,maxbond,nloc,udat,wfnfile,maxcritov,
     &      critov,ixcent,largwr,warn,skiph,stdout,stderr)
          else
            call ulewisoqs (aomnat,mal,mbe,nfrag,nwithin,ifrag,
     &      mxcent,covx,epsneg,maxbond,nloc,udat,wfnfile,maxcritov,
     &      critov,ixcent,largwr,warn,skiph,stdout,stderr)
          endif
        endif
        call timestamp (line)
        write (stdout,2000) line(1:leng(line))
        stop 'edf: End of Run Ok !!!'
      else
        if (onlynat) then
          write (stderr,322)
          stop
        endif
      endif
c
c.....Analyze which of the canonical MOs can be considered as 
c     CORE orbitals, and determine the fragment they belong to.
c
      if (.not.allocated(mocogrp)) then
        allocate (mocogrp(ngroup),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate mocogrp()'
      endif
      if (.not.allocated(mocogrpa)) then
        allocate (mocogrpa(ngroup),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate mocogrpa()'
      endif
      if (.not.allocated(mocogrpb)) then
        allocate (mocogrpb(ngroup),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate mocogrpb()'
      endif
      if (.not.allocated(icogrp)) then
        allocate (icogrp(ngroup,nmo),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate icogrp()'
      endif
      mocogrp  = izero
      mocogrpa = izero
      mocogrpb = izero
      mocore  = izero
      if (corecan) then
        call coremos (sg,ngroup,nmo,ifilc,
     &  ndets,ncore,nelcorea,nelcoreb,nact,nelact,nel,nalpha,nbeta,
     &  ovcritcan,mocore,mocogrp,mocogrpa,mocogrpb,icore,icorea,
     &  icoreb,moval,ival,stdout,ok,largwr,icorr,largeline,wfnfile)
c
c-------write NEW AOM and WFN files with the core MOs elliminated
c
        oldlarge = largeline
        ln = leng(largeline)
        if (ok) then
          largeline(1:4) = 'AOM '
          call suppress 
     &      (aom,udat,filedat,wfnfile,largeline(1:ln),stdout)
          largeline = oldlarge
          largeline(1:4) = 'WFN '
          call suppress 
     &      (aom,udat,filedat,wfnfile,largeline(1:ln),stdout)
        else
          write (stdout,641)
        endif
      endif
c
c.....Modify the minimum and maximum values of ALPHA and BETA electrons 
c     discounting the core orbitals.
c
      if (.not.(skipcan.and.skiploc.and.skiplocx)) write (stdout,504) 
      nfutgrp=0
      do i=1,ngroup
         eleca(2,i)=eleca(2,i)-nelcorea
         elecb(2,i)=elecb(2,i)-nelcoreb
         if (.not.(skipcan.and.skiploc.and.skiplocx)) then
           write (stdout,500) i,nfugrp(i),(ifugrp(j,i),j=1,nfugrp(i))
           write (stdout,502) ' # ALPHA ',eleca(1,i),eleca(2,i)
           write (stdout,502) ' # BETA  ',elecb(1,i),elecb(2,i)
         endif
         nfutgrp=nfutgrp+nfugrp(i)
      enddo
      write (stdout,*) 
      if (nfutgrp.ne.ncent) then
        write (stderr,*)
        write (stderr,*) '# Basins not present in any group ?'
        write (stderr,*)
        stop
      endif

      naval = mal-nelcorea
      nbval = mbe-nelcoreb
      nval  = nel-nelcorea-nelcoreb

      if (nindep.and.ngroup.gt.2) then
        call configs (numdet,nmo,nelact,nact,ifilc,stderr,wfnfile)
        call indepgr (epsdet,probcut,numdet,nelact,nact,ncent,nmo,ncore,
     &  ngroup,mal,mbe,mocore,nel,moval,ifilc,stdout,stderr,aom,nfugrp,
     &  ifugrp,mocogrp,ival,doentropy,orderp,short,wfnfile)
      endif
      if (bonding) then
        call configs (numdet,nmo,nelact,nact,ifilc,stderr,wfnfile)
        call optedf (epsdet,probcut,numdet,nelact,nact,ncent,nmo,ncore,
     &       ngroup,mal,mbe,mocore,nel,moval,ifilc,stdout,stdin,aom,sg,
     &       nfugrp,ifugrp,mocogrp,ival,short,wfnfile)
      endif
      if (maxwfn) then
        call configs (numdet,nmo,nelact,nact,ifilc,stderr,wfnfile)
        call wfnmax (stdin,stdout,stderr)
      endif
c
c.....Exact computation of probabilities.
c
      if (.not.skipcan) then
c
c.......Obtain the resonance structures
c
        call timer (2,ipid1,'_rnprobs  ',-1)
c
c-------Determine the number of ALPHA and BETA probabilities
c
c-------AQUI hay un error: Cuando hay orbitales de CORE con occupcion
c       1.0 ===> nproba y nprobb se calculan MAL
c
        nproba=1
        nprobb=1
        mmm=1
        malco=mal-nelcorea
        mbeco=mbe-nelcoreb
        do i=1,ngroup-1
          nproba=nproba*(malco+i)
          nprobb=nprobb*(mbeco+i)
          mmm=mmm*i
        enddo
        nproba=nproba/mmm
        nprobb=nprobb/mmm
c
c.......Allocate the arrays of alpha and beta RSRSs.
c  
        if (.not.allocated(resnca)) then
          allocate (resnca(nproba,ngroup),stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot allocate resnca()'
        endif
        if (.not.allocated(resncb)) then
          allocate (resncb(nprobb,ngroup),stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot allocate resncb()'
        endif
c
        
        call rnprobs (eleca,resnca,malco,ngroup,nproba,stdout)
        call rnprobs (elecb,resncb,mbeco,ngroup,nprobb,stdout)
        call timer (4,ipid1,'_rnprobs  ',-1)

        write (stdout,1130) randommin,randommax
        if (nalpha.eq.nbeta) then
          if (ngroup.eq.2.and.recur.and.(mal.eq.mbe)) then
            call rcalcedf  (epsdet,probcut,nproba,nprobb,nprob,ngroup,
     &        numdet,nmo,ncore,nelact,nact,nel,moval,ival,mocogrp,
     &        mocogrpa,mocogrpb,naval,ifilc,stdout,stderr,sg,
     &        resnca,resncb,mulliken,largwr,doentropy,orderp)
          else
            call binedf (randommin,randommax,
     &        epsdet,probcut,epsproba,nproba,nprobb,
     &        nprob,ngroup,numdet,nmo,ncore,nelact,nel,moval,ival,
     &        mocogrp,mocogrpa,mocogrpb,
     &        naval,nbval,ifilc,stdout,stderr,wfnfile,sg,
     &        resnca,resncb,largwr,doentropy,orderp,memedf)
          endif
        else
         call binedf (randommin,randommax,
     &      epsdet,probcut,epsproba,nproba,nprobb,
     &      nprob,ngroup,numdet,nmo,ncore,nelact,nel,moval,ival,
     &      mocogrp,mocogrpa,mocogrpb,
     &      naval,nbval,ifilc,stdout,stderr,wfnfile,sg,
     &      resnca,resncb,largwr,doentropy,orderp,memedf)
        endif
        if (allocated(resnca)) then
          deallocate (resnca,stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot deallocate resnca()'
        endif
        if (allocated(resncb)) then
          deallocate (resncb,stat=ier)
          if (ier.ne.0) stop 'edf.f: Cannot deallocate resncb()'
        endif
*       call timer (5,ipid,'_edf      ',stdout)
*       call timestamp (line)
*       write (stdout,2000) line(1:leng(line))
*       stop 'edf: END OF RUN Ok !!!'
      endif
c
c-----Fragment Natural Orbitals analysis
c
      if (dofno) then
         do i=1,nfno
*          call fno (aom,inside(1,i),rhf,rohf,uhf,cciqa,icorr,mal,mbe,
*    &      ninside(i),nprims,nmo,ncent,epseigen,udat,wfnfile,stdout,
*    &      stderr,largwr)

*          call otherfno (aom,inside(1,i),rhf,rohf,uhf,cciqa,icorr,mal,
*    &        mbe,ninside(i),nprims,nmo,ncent,epseigen,udat,wfnfile,
*    &        stdout,stderr,largwr)

*          call newfno (aom,inside(1,i),rhf,rohf,uhf,cciqa,icorr,mal,
*    &        mbe,ninside(i),nprims,nmo,ncent,epseigen,udat,wfnfile,
*    &        stdout,stderr,largwr)

           call fnofrag (aom,inside(1,i),cciqa,mal,mbe,ninside(i),
     &        nprims,nmo,ncent,udat,wfnfile,stdout,stderr,largwr)

*          call fnoeos (aom,nfugrp,ifugrp,rhf,rohf,uhf,cciqa,icorr,
*    &        mal,mbe,nprims,nmo,ncent,ngroup,epseigen,udat,wfnfile,
*    &        stdout,stderr,largwr)

*          call nnn (aom,inside(1,i),rhf,rohf,uhf,cciqa,icorr,mal,
*    &        mbe,ninside(i),nprims,nmo,ncent,epseigen,udat,wfnfile,
*    &        stdout,stderr,largwr)

         enddo
         deallocate (inside)
         deallocate (ninside)
         call timer (5,ipid,'_edf      ',stdout)
         call timestamp (line)
         write (stdout,2000) line(1:leng(line))
         stop 'edf: End of Run Ok !!!'
      endif
c
c.....Localization or delocalization is performed
c
      if (.not.allocated(vrot)) then
        allocate (vrot(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate vrot()'
      endif
      if (.not.allocated(canmo)) then
        allocate (canmo(nmo,nmo),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate canmo()'
      endif
c
c-----Localization of the Canonical MOs
c
      lewis = canloc
      lewis = lewis  .or. chemiloc
      lewis = lewis  .or. dafhloc
      lewis = lewis  .or. aomloc
      lewis = lewis  .or. otherdafh
      lewis = lewis  .or. openloc
      lewis = lewis  .or. otheropen
      if (lewis) then
        if (.not.canloc) allcloc=.true.
        if (allcloc) then
          if (.not.allocated(cc)) then
            allocate (cc(nmo,nprims),stat=ier)
            if (ier.ne.0) stop 'edf.f: Cannot allocate cc()'
          endif
c
c---------Normalize MOs
c
          cc=coef(1+nmo:nmo+nmo,:)
          call normalizemos (cc,nmo,nprims)
c
c---------Multipole moments of canonical MOs are obtained.
c
          write (stdout,334) 'Canonical MOs'
          call newdrvmult (nmo,nprims,ncent,cc,stdout)
          sgallc(1:ncent,1:nmo,1:nmo)=aom(1:ncent,1:nmo,1:nmo)
          if (mod(nel,2).eq.0) then
            nelhalf=nel/2
          else
            nelhalf=(nel+1)/2
          endif
          if (chemiloc) then
             woccup = .false.
             call uhfchemloc (sgallc,maxcritov,critov,damps,critics,
     &         covx,hesel,stdout,largwr,.true.,wfnfile)
               goto 1234
          elseif (aomloc) then
             woccup = .false.
             call uhfitaomloc (sgallc,mxcent,maxcritov,critov,largwr,
     &         lw,wfnfile)
               goto 1234
          elseif (dafhloc.or.otherdafh) then
c
c...........In case that the WFN is a multideterminantal one, we need to 
c           compute the DAFH in every atom by integrating rho^XC.
c           
            if (icorr.or.cciqa) then
              if (nalpha.eq.nbeta) then
                allocate (dafh(ncent,nmo,nmo))
                call compdafh (aom,dafh,nmo,ncent,stderr,wfnfile)
              else
                stop '# Unable DAHFLOC with multidet open shell WFNs'
              endif
              woccup = .false.
              if (dafhloc) call itdafhc (sgallc,dafh,vrot,mxcent,covx,
     &             maxbond,nloc,nelhalf,udat,wfnfile,maxcritov,
     &             critov,ixcent,largwr,warn,skiph,stdout)
              if (otherdafh) call dafhdo (sgallc,dafh,vrot,nloc,
     &            nelhalf,udat,wfnfile,largwr,warn,lw)
                deallocate (dafh)
            else
              if (dafhloc) call itdafh (sgallc,vrot,mxcent,covx,
     &          maxbond,udat,ifilc,wfnfile,maxcritov,critov,dafhselect,
     &          ndafhsel,idafhsel,ixcent,largwr,warn,skiph,stdout)
              if (otherdafh) call dafhdo (sgallc,sgallc,vrot,nloc,
     &            nelhalf,udat,wfnfile,largwr,warn,lw)
              nloc=nmo
              if (warn) goto 3333
            endif
          elseif (openloc) then
            if (uhf) then
              call uhfdrvloc (sgallc,mxcent,covx,epsneg,maxbond,mal,
     &          mbe,nloc,udat,wfnfile,maxcritov,critov,ixcent,largwr,
     &          warn,skiph,stdout,stderr)
            endif
            if (rhf.or.(icorr.and.mal.eq.mbe)) then
              call casoqsloc (sgallc,mxcent,covx,epsneg,maxbond,
     &          nloc,mal,udat,wfnfile,maxcritov,critov,
     &          ixcent,largwr,warn,skiph,stdout,stderr)
              if (nloc.ne.mal) then
                stop '# edf.f: Failed localizing Canonical MOs'
              endif
            endif
            if (rohf.or.(icorr.and.mal.ne.mbe)) then
              call casuhfdrvloc (sgallc,mxcent,covx,epsneg,
     &          maxbond,nloc,mal,mbe,udat,wfnfile,maxcritov,
     &          critov,ixcent,largwr,warn,skiph,stdout,stderr)
            endif
            if (cciqa.and.mal.eq.mbe) then
              call casoqsloc (sgallc,mxcent,covx,epsneg,maxbond,
     &          nloc,mal,udat,wfnfile,maxcritov,critov,
     &          ixcent,largwr,warn,skiph,stdout,stderr)
              if (nloc.ne.mal) then
                stop '# edf.f: Failed localizing Canonical MOs'
              endif
            endif
            if (cciqa.and.mal.ne.mbe) then
              stop '# edf.f: Open-shell CC WFNs not allowed'
            endif
            goto 1234
          elseif (otheropen) then
            if (uhf) then
              call udrvother (sgallc,mxcent,covx,epsneg,maxbond,
     &          nloc,mal,mbe,udat,wfnfile,maxcritov,critov,ixcent,
     &          largwr,warn,skiph,stdout,stderr)
            endif
            if (rhf.or.(icorr.and.mal.eq.mbe)) then
              call drvother (sgallc,mxcent,covx,epsneg,maxbond,
     &          nloc,mal,mbe,udat,wfnfile,maxcritov,critov,ixcent,
     &          largwr,warn,skiph,stdout,stderr)
            endif
            if (rohf.or.(icorr.and.mal.ne.mbe)) then
              call udrvother (sgallc,mxcent,covx,epsneg,maxbond,
     &          nloc,mal,mbe,udat,wfnfile,maxcritov,critov,ixcent,
     &          largwr,warn,skiph,stdout,stderr)
            endif
            if (cciqa.and.mal.eq.mbe) then
              call drvother (sgallc,mxcent,covx,epsneg,maxbond,
     &          nloc,mal,mbe,udat,wfnfile,maxcritov,critov,ixcent,
     &          largwr,warn,skiph,stdout,stderr)
            endif
            if (cciqa.and.mal.ne.mbe) then
              stop '# edf.f: Open-shell CC WFNs not allowed'
            endif
            goto 1234
          else
            call canlmo (sgallc,vrot,nmo,ncent,critics,
     &            stdout,largwr,.true.)
          endif
c
c---------Test consistency of the obtained AOM elements.
c
          call aomtest 
     &       (sgallc(:,1:nloc,1:nloc),nloc,ncent,critics,stdout,largwr)
c
          do i=1,nloc
            do ip=1,nprims
              cc(i,ip) = dot_product(coef(nmo+1:nmo+nmo,ip),vrot(i,:))
            enddo
          enddo
c
c---------AOM matrices are computed with the localized MOs
c
          coef(1:nloc,1:nprims)=cc(1:nloc,1:nprims)
          if (.not.(icorr.or.cciqa)) then
            if (mulliken) then
              call aomulliken (sgallc,1,.false.,stdout,wfnfile)
            elseif (mindef) then
              call aomindef (sgallc,1,.false.,stdout,wfnfile)
            elseif (mindefrho.or.netrho.or.promrho.or.hesel) then
              call aomindefrho (sgallc,alphahesel,betahesel,
     &         irho,nrad,nang,irmesh,rhopow,.false.,coef,stdout,wfnfile)
            elseif (lowdin) then
              call aomlowdin (sgallc,1,.false.,stdout,wfnfile)
            elseif (becke) then
              call aombecke (sgallc,nrad,nang,pow,irmesh,.true.,
     &        cc,stdout,wfnfile)
            endif
          endif
c
c---------Multipole moments of localized MOs are obtained.
c
          write (stdout,334) 'Localized MOs'
          call newdrvmult (nloc,nprims,ncent,cc(1:nloc,:),stdout)
          inmo=nloc+nloc
          inpr=nprims
          if (chemiloc) then
            wfnloc=trim(wfnfile)//"chem"
          elseif (aomloc) then
            wfnloc=trim(wfnfile)//"aom"
          elseif (dafhloc) then
            wfnloc=trim(wfnfile)//"dafh"
          else
            wfnloc=trim(wfnfile)//"loc"
          endif
          open (unit=udatnw,file=wfnloc,status='unknown')
          write (stdout,351) wfnloc(1:leng(wfnloc))
c
          if (icorr.or.cciqa) then
            allocate (oc(nloc),stat=ier)
            if (ier.ne.0) stop 'edf.f: Cannot allocate oc()'
            oc(1:nloc)=occ(1:nloc)
            call cdafhmos 
     &           (udat,udatnw,stderr,oc,cc(1:nloc,:),nloc,inpr)
            deallocate (oc,stat=ier)
            if (ier.ne.0) stop 'edf.f: Cannot deallocate oc()'
          else
            woccup = .false.
            call wrtwfn (udat,udatnw,stderr,ifilc,coef,inmo,inpr)
c
c-----------Fill in AOM file with localized MOs when any of dafhloc,
c           otherdafh, openloc, or otheropen are .true.
c
            if (dafhloc.or.otherdafh) then
              lu19    =  19
              exfil   = .true.
              fileloc = trim(wfnloc)//'.aom'
              do while (exfil)
                inquire (unit=lu19,opened=exfil)
                if (.not.exfil) then
                  open (unit=lu19,file=trim(fileloc),iostat=ios)
                  if (ios.ne.0) then
                    write (0,*) 
     &              ' # edf.f: Error openning '//trim(fileloc)
                    stop
                  endif
                else
                  lu19=lu19+1
                endif
              enddo
              if (dafhloc.or.otherdafh) then
                write (stdout,3511) trim(fileloc),'DAF'
              else
                write (stdout,3511) trim(fileloc),'OPENLOC'
              endif
              do i=1,ncent
                write (lu19,'(I4,a)') i,' <--- AOMLOC in this center'
                write (lu19,80) ((sgallc(i,m,j),m=1,j),j=1,nmo)
              enddo
            endif
          endif
          close (udatnw)
c
 3333     continue
          if (allocated(cc)) then
            deallocate (cc,stat=ier)
            if (ier.ne.0) stop 'edf.f: Cannot deallocate cc()'
          endif
        else
          sgloc=zero
          do ngr=1,ngroup
            do i=1,nfugrp(ngr)
              sgloc(ngr,:,:) = sgloc(ngr,:,:) + aom(ifugrp(i,ngr),:,:)
            enddo
          enddo
          call canlmo 
     &         (sgloc,vrot,nmo,ngroup,critics,stdout,largwr,.true.)
        endif
      endif
 1234 continue
c
      if (.not.allocated(oc)) then
        allocate (oc(nmo),stat=ier)
        if (ier.ne.0) stop 'edf.f: Cannot allocate oc()'
      endif
      oc(1:nmo)=occ(1:nmo)
c
c.....If the isopycnic localization is performed on all the centers,
c     array sgallc() represents the AOM between localized MOs. Before
c     calling the localization routine (ruedmis), it is initialized to
c     the AOM between the natural MOs,
c
      if (localize) then
         if (allcloc) then
           allocate (sgtmp(nmo,nmo))
           allocate (vtmp(nmo))
           sgallc=zero
           do igr=1,ncent
             sgtmp=aom(igr,:,:)
             do k=1,nmo
               vtmp(:)=matmul(sgtmp,v1mata(:,k))
               do i=1,nmo
                 sgallc(igr,i,k)=dot_product(v1mata(:,i),vtmp(:))
               enddo
             enddo
           enddo
           deallocate (sgtmp)
           deallocate (vtmp)
         else
           sgloc(1:ngroup,1:nmo,1:nmo)=sgnat(1:ngroup,1:nmo,1:nmo)
         endif
         if (.not.allocated(cc)) then
           allocate (cc(nmo+nmo,nprims),stat=ier)
           if (ier.ne.0) stop 'edf.f: Cannot allocate cc()'
         endif
         cc(nmo+1:nmo+nmo,1:nprims)=coef(nmo+1:nmo+nmo,1:nprims)
c
c--------Measuring the locality of Localized MOs. First, Normalize MOs
c
         call normalizemos (cc(1+nmo:nmo+nmo,:),nmo,nprims)
         write (stdout,334) 'Canonical MOs'
         call newdrvmult (nmo,nprims,ncent,cc(1+nmo:nmo+nmo,:),stdout)
c
c........Perform an isopycnic localization.
c
         if (allcloc) then    ! In all centers
            write (stdout,1121)
            vrot=0d0
            call ruedmis (sgallc,vrot,oc,nmo,ncent,stdout,largwr,.true.)
c
c...........Obtain group overlap integrals sgloc() between localized MOs 
c           from sgallc(), the AOM between these orbitals
c
            sgloc=zero
            if (ngroup.eq.1) then
              forall (m=1:nmo) sgloc(1,m,m)=one
            else
              do k=1,ncent
                do ngr=1,ngroup
                  do i=1,nfugrp(ngr)
                    if (ifugrp(i,ngr).eq.k) then
                       sgloc(ngr,:,:)=sgloc(ngr,:,:)+sgallc(k,:,:)
                    endif
                  enddo
                enddo
              enddo
            endif
c
c-----------Fill in AOM file with localized MOs
c
            lu19    =  19
            exfil   = .true.
            fileloc = trim(wfnfile)//'loc.aom'
            do while (exfil)
              inquire (unit=lu19,opened=exfil)
              if (.not.exfil) then
                open (unit=lu19,file=trim(fileloc),iostat=ios)
                if (ios.ne.0) then
                  write (0,*) ' # edf.f: Error openning '//trim(fileloc)
                  stop
                endif
              else
                lu19=lu19+1
              endif
            enddo
            write (stdout,3510) trim(fileloc)
            if (ioverlap.eq.2.or.ioverlap.eq.3) then
              write (lu19,'(3I5)') nbasins, nmo, natoms
              write (lu19,'(1000(1x,a))') (trim(namebas(i)),i=1,ncent)
              do i=1,ncent
                write (lu19,80) ((sgallc(i,j,m),m=j,nmo),j=1,nmo)
              enddo
            else
              do i=1,ncent
                write (lu19,'(I4,a)') i,' <--- AOMLOC in this center'
                write (lu19,80) ((sgallc(i,m,j),m=1,j),j=1,nmo)
              enddo
            endif
         else                      ! In all groups
            write (stdout,112)
            call ruedmis (sgloc,vrot,oc,nmo,ngroup,stdout,largwr,.true.)
         endif
c
c........Expand the localized MOs (LMO) in terms of primitive Gaussians 
c        ( cc() array), and obtain the transformation matrix from 
c        canonical orbitals to localized orbitals. This matrix is stored
c        in the canmo() array.
c        
         cc(nmo+1:nmo+nmo,1:nprims)=coef(nmo+1:nmo+nmo,1:nprims)
         do i=1,nmo
           do j=1,nmo
             temp=zero
             do k=1,nmo
               temp=temp+v1mata(j,k)*vrot(i,k)
             enddo
             canmo(j,i)=temp
           enddo
           do j=1,nprims
             temp = 0d0
             do k=1,nmo
               temp = temp + canmo(k,i) * coef(k+nmo,j)
             enddo
             cc(i,j) = temp
           enddo
         enddo
c
c--------Measuring the locality of Localized MOs. First, Normalize MOs
c
         call normalizemos (cc(1:nmo,:),nmo,nprims)
         write (stdout,334) 'Localized MOs'
         call newdrvmult (nmo,nprims,ncent,cc(1:nmo,:),stdout)
c
c........Obtain the transformation matrix from LMOs to CANMOs, 
c        i.e. the inverse of canmo(k,i).
c
         if (.not.allocated(work)) then
           allocate (work(nmo),stat=ier)
           if (ier.ne.0) stop 'edf.f: Cannot allocate work()'
         endif
         if (.not.allocated(ipvt)) then
           allocate (ipvt(nmo),stat=ier)
           if (ier.ne.0) stop 'edf.f: Cannot allocate ipvt()'
         endif
c
         call dgeco (canmo,nmo,nmo,ipvt,rcond,work)
         job=01
         deter=zero
         call dgedi (canmo,nmo,nmo,ipvt,deter,work,job)
c
         if (allocated(work)) then
           deallocate (work,stat=ier)
           if (ier.ne.0) stop 'edf.f: Cannot deallocate work()'
         endif
c
         if (allocated(ipvt)) then
           deallocate (ipvt,stat=ier)
           if (ier.ne.0) stop 'edf.f: Cannot deallocate ipvt()'
         endif
c
         if (largwr) then
           write (stdout,41)
           do m=1,nmo
             write (stdout,331) m
             write (stdout,30) (canmo(j,m),j=1,nmo)
           enddo
         endif
c
c........lmowfn.f routine requires the transpose of canmo().
c
         do i=2,nmo
           do j=1,i-1
             temp=canmo(i,j)
             canmo(i,j)=canmo(j,i)
             canmo(j,i)=temp
           enddo
         enddo
c
c........Write a WFN file changing the input natural MOs by the recently
c        computed localized MOs.
c
         inmo=nmo
         inpr=nprims
         occ(1:nmo)=oc(1:nmo)
         wfnloc=wfnfile(1:leng(wfnfile))//"loc"
         open (unit=udatnw,file=wfnloc,status='unknown')
         write (stdout,351) wfnloc(1:leng(wfnloc))
         woccup = .true.
         call wrtwfn (udat,udatnw,stderr,ifilc,cc(1:nmo,:),inmo,inpr)
         if (allocated(cc)) then
           deallocate (cc,stat=ier)
           if (ier.ne.0) stop 'edf.f: Cannot deallocate cc()'
         endif
         close (udatnw)
         write (stdout,'(1x,a)') '# Done'
c
         if ((.not.skiploc).or.(.not.skiplocx)) then
           if (.not.allocated(ivaloc)) then
             allocate (ivaloc(nmo),stat=ier)
             if (ier.ne.0) stop 'edf.f: Cannot allocate ivaloc()'
           endif
c
c..........Analyze which of the localized MOs may be considered as CORE
c          MOs, and determine to which fragment they belong.
c
c AQUI
           mocogrp=izero
           if (coreloc) then
             do i=1,ngroup
               eleca(2,i)=eleca(2,i)+mocore
               elecb(2,i)=elecb(2,i)+mocore
             enddo
             mocore=izero
             do ig=1,ngroup
               do i=1,nmo
                 if (sgloc(ig,i,i).gt.ovcritloc) then
                   mocore=mocore+1
                   mocogrp(ig)=mocogrp(ig)+1
                   icore(mocore)=i
                   icogrp(ig,mocogrp(ig))=i
                 endif
               enddo
             enddo
             movaloc=izero
             do i=1,nmo
               j=1
               iscore=.false.
               do while (j.le.mocore.and.(.not.iscore))
                 if (i.eq.icore(j)) iscore=.true.
                 j=j+1
               enddo
               if (.not.iscore) then
                 movaloc=movaloc+1
                 ivaloc(movaloc)=i
               endif
             enddo
             write (stdout,223) movaloc,(ivaloc(i),i=1,movaloc)
             write (stdout,221) ovcritloc
             do ig=1,ngroup
               mocg=mocogrp(ig)
               write (stdout,222) ig,mocg,(icogrp(ig,k),k=1,mocg)
               do k=1,mocg
                 write (stdout,225) icogrp(ig,k),icogrp(ig,k),ig,
     &           sgloc(ig,icogrp(ig,k),icogrp(ig,k))
               enddo
             enddo
c
             nval=nel-2*mocore
             naval=mal-mocore
             nbval=mbe-mocore
c
             do i=1,ngroup
               eleca(2,i)=eleca(2,i)-mocore
               elecb(2,i)=elecb(2,i)-mocore
             enddo
c
             call timer (2,ipid1,'_rnprobs  ',-1)
c
c............Allocate the arrays resnca() and resncb() computing firt
c            the number of resonance structures.
c
             nproba=1
             nprobb=1
             mmm=1
             do i=1,ngroup-1
               nproba=nproba*(naval+i)
               nprobb=nprobb*(nbval+i)
               mmm=mmm*i
             enddo
             nproba=nproba/mmm
             nprobb=nprobb/mmm
             if (.not.allocated(resnca)) then
               allocate (resnca(nproba,ngroup),stat=ier)
               if (ier.ne.0) stop 'edf.f: Cannot allocate resnca()'
             endif
             if (.not.allocated(resncb)) then
               allocate (resncb(nprobb,ngroup),stat=ier)
               if (ier.ne.0) stop 'edf.f: Cannot allocate resncb()'
             endif
c
             call rnprobs (eleca,resnca,naval,ngroup,nproba,stdout)
             call rnprobs (elecb,resncb,nbval,ngroup,nprobb,stdout)
             call timer (4,ipid1,'_rnprobs  ',-1)
           else
             mocore=izero
             movaloc=nmo
             nval=nel
             naval=mal
             nbval=mbe
             forall (i=1:movaloc) ivaloc(i)=i
             do i=1,ngroup
               eleca(2,i)=eleca(2,i)+mocore
               elecb(2,i)=elecb(2,i)+mocore
             enddo
             call timer (2,ipid1,'_rnprobs  ',-1)
c
c............Allocate the arrays resnca() and resncb() computing firt
c            the number of resonance structures.
c
             nproba=1
             nprobb=1
             mmm=1
             do i=1,ngroup-1
               nproba=nproba*(naval+i)
               nprobb=nprobb*(nbval+i)
               mmm=mmm*i
             enddo
             nproba=nproba/mmm
             nprobb=nprobb/mmm
             if (.not.allocated(resnca)) then
               allocate (resnca(nproba,ngroup),stat=ier)
               if (ier.ne.0) stop 'edf.f: Cannot allocate resnca()'
             endif
             if (.not.allocated(resncb)) then
               allocate (resncb(nprobb,ngroup),stat=ier)
               if (ier.ne.0) stop 'edf.f: Cannot allocate resncb()'
             endif
c
             call rnprobs (eleca,resnca,naval,ngroup,nproba,stdout)
             call rnprobs (elecb,resncb,nbval,ngroup,nprobb,stdout)
             call timer (4,ipid1,'_rnprobs  ',-1)
           endif
c
c..........Find expansion coefficients of the wave function in terms
c          of Slater dets built in with localized MOs.
c
           call lmowfn (canmo,epswfnloc,numdet,nmo,ncore,nelact,
     &        nel,mal,mbe,ndetlmo,largwr,wfnfile,cicoef,ifilc,ifilmo,
     &        stdout,stderr)
         endif
c
c........Exact calculation of probabilities with Localized MOs
c
         if (.not.skiploc) then
           write (stdout,1131) randommin,randommax
           if (numdet.eq.1) then
             call nbasab (probcut,nproba,nprobb,nprob,ngroup,nmo,
     &         movaloc,ivaloc,mocogrp,nel,nval,naval,nbval,ifilc,
     &         stdout,sgloc,resnca,resncb,largwr,doentropy,orderp)
           else
             call binedf (randommin,randommax,
     &         epsdet,probcut,epsproba,nproba,nprobb,
     &         nprob,ngroup,ndetlmo,nmo,0,nel,nel,movaloc,ivaloc,
     &         mocogrp,mocogrpa,mocogrpb,
     &         naval,nbval,ifilmo,stdout,stderr,wfnfile,sgloc,
     &         resnca,resncb,largwr,doentropy,orderp,memedf)
           endif
         endif       
c
c........Approximate calculation of probabilities with Localized MOs
c
         if (.not.skiplocx) then
            if (.not.allocated(resnca)) then
              allocate (resnca(nproba,ngroup),stat=ier)
              if (ier.ne.0) stop 'edf.f: Cannot allocate resnca()'
            endif
            if (.not.allocated(resncb)) then
              allocate (resncb(nprobb,ngroup),stat=ier)
              if (ier.ne.0) stop 'edf.f: Cannot allocate resncb()'
            endif
           if (.not.icorr) then
             call xnbasab (probcut,nproba,nprobb,nprob,ngroup,nmo,
     &            movaloc,ivaloc,mocogrp,nel,nval,naval,nbval,ifilc,
     &            stdout,sgloc,resnca,resncb,largwr,doentropy,orderp)
           else
             if (nalpha.eq.nbeta) then
               call xcalcedf  (epsdet,probcut,nproba,nprobb,nprob,
     &              ngroup,ndetlmo,wfnfile,nmo,0,nel,nel,
     &              movaloc,ivaloc,mocogrp,naval,ifilmo,stdout,
     &              stderr,sgloc,resnca,resncb,largwr,doentropy,orderp)
             else
               call xcalcedfd (epsdet,probcut,nproba,nprobb,nprob,
     &              ngroup,ndetlmo,wfnfile,nmo,0,nel,nel,
     &              movaloc,ivaloc,mocogrp,naval,nbval,ifilmo,stdout,
     &              sgloc,resnca,resncb,largwr,doentropy,orderp)
             endif
           endif
         endif
c
         if (allocated(ivaloc)) then
           deallocate (ivaloc,stat=ier)
           if (ier.ne.0) then
             stop 'edf.f: Cannot deallocate ivaloc()'
           endif
         endif
         if (allocated(resnca)) then
           deallocate (resnca,stat=ier)
           if (ier.ne.0) stop 'edf.f: Cannot deallocate resnca()'
         endif
         if (allocated(resncb)) then
           deallocate (resncb,stat=ier)
           if (ier.ne.0) stop 'edf.f: Cannot deallocate resncb()'
         endif
      endif
c
c-----Probabilities for specific resonance structures.
c
      if (nrsrs.gt.0) then
        if (ngroup.eq.0) then
          stop ' # edf.f: Number of groups not read in yet'
        endif
        call configs (numdet,nmo,nelact,nact,ifilc,stderr,wfnfile)
        do i=1,nrsrs
          call probres (epsdet,numdet,ok,ngroup,mal,mbe,nmo,
     &                ncore,nel,stdout,stderr,maxalpha,minalpha,
     &                maxbeta,minbeta,probal,nres(i))
          if (ok) write (stdout,116) probal,nres(i)(1:leng(nres(i)))
        enddo
      endif
c
c.....All two-center DIs are computed.
c
      if (ditwocent) then
        call configs (numdet,nmo,nelact,nact,ifilc,stderr,wfnfile)
        call twocendi (aom,epsbond,epsdet,probcut,numdet,ncent,nmo,
     &   nel,nval,ncore,nelact,moval,mal,mbe,ifilc,stdout,stderr,ngroup,
     &   ival,wfnfile,largwr)
      endif
c
      close (ifilc,status='delete')
      close (ifilmo,status='delete')
c
      call timer (5,ipid,'_edf      ',stdout)
      call timestamp (line)
      write (stdout,2000) line(1:leng(line))
      write (stderr,'(a)') 
     & " # End of Run Ok for WFN file "//trim(wfnfile)
c
c.....Formats
c
 1000 format(/,1x,101('+'),/,
     &  1x,'|',3('  ______    ______    ______   '),6x,'|',/,
     &  1x,'|',3(' /_____/\  /_____/\  /_____/\  '),6x,'|',/,
     &  1x,'|',3(' \::::_\/_ \:::_ \ \ \::::_\/_ '),6x,'|',/,
     &  1x,'|',3('  \:\/___/\ \:\ \ \ \ \:\/___/\'),6x,'|',/,
     &  1x,'|',3('   \::___\/  \:\ \ \ \ \:::._\/'),6x,'|',/,
     &  1x,'|',3('    \:\____/\ \:\/.:| | \:\ \  '),6x,'|',/,
     &  1x,'|',3('     \_____\/  \____/_/  \_\/  '),6x,'|',/,
     &  1x,'|',3('                               '),6x,'|',/,
     &  1x,101('+'),/,31x,'ELECTRON NUMBER DISTRIBUTION FUNCTIONS',/,
     & 31x,'(c) E. Francisco & A. Martín Pendas, 2022',/,37x,
     &  'University of Oviedo',/,1x,101('+'),/,21x,
     & 'Calculation starts on ',a,/,1x,101('+'))
*1000 format(/,
*    &' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/
*    &'  ______    ______    ______    ______    ______    ______   ',/
*    &' /_____/\  /_____/\  /_____/\  /_____/\  /_____/\  /_____/\  ',/
*    &' \::::_\/_ \:::_ \ \ \::::_\/_ \::::_\/_ \:::_ \ \ \::::_\/_ ',/
*    &'  \:\/___/\ \:\ \ \ \ \:\/___/\ \:\/___/\ \:\ \ \ \ \:\/___/\',/
*    &'   \::___\/  \:\ \ \ \ \:::._\/  \::___\/_ \:\ \ \ \ \:::._\/',/
*    &'    \:\____/\ \:\/.:| | \:\ \     \:\____/\ \:\/.:| | \:\ \  ',/
*    &'     \_____\/  \____/_/  \_\/      \_____\/  \____/_/  \_\/  ',/
*    &'                                                             ',/
*    &' ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/
*    &  '                ELECTRON NUMBER DISTRIBUTION FUNCTIONS     ',/,
*    &  '              (c) E. Francisco & A. Martín Pendas, 2022    ',/,
*    &  '                      University of Oviedo                 ',/,
*    &  ' Calculation starts on ',a,/,1x,60('+'))
 2000 format (' # ',54('+'),/,
     &        ' # Calculation ends on ',a,/,
     &        ' # ',54('+'))
 321  format (' # edf.f: Bad format in record #1 of input file',/,
     .        ' # edf.f: IOVERLAP variable not read in')
 33   format (1x,'# AOM file                          = ',a)
 330  format (1x,'# CRITIC2 AOM file = ',a)
 34   format (1x,'# Wave Function file                = ',a)
 36   format (1x,'# Number of determinants in WFN     = ',I10)
 363  format (1x,'# Type of single-determinant WFN    = ',7x,a)
 361  format (1x,'# Number of determinants in WFN     = ',I10,/,
     &        1x,'# Number of active electrons in WFN = ',I10,/,
     &        1x,'# Number of core orbitals in WFN    = ',I10,/,
     &        1x,'# Number of active orbitals in WFN  = ',I10)
 37   format (1x,'# EPSDET                            = ',E18.10,/,
     &        1x,'# EPSWFN                            = ',E18.10,/,
     &        1x,'# EPSWFNLOC                         = ',E18.10,/,
     &        1x,'# PROBCUT                           = ',E18.10,/,
     &        1x,'# EPSPRBA                           = ',E18.10)
 55   format (1x,'# AOM read in with ',a,' format')
 56   format (1x,'# Analysis ',a,' including SPIN probabilities')
 57   format (1x,'# Complex Atomic Overlap Matrix is read in')
 58   format (1x,'# Complex Group Overlap integrals are computed')
 103  format (1x,A4,4x,I3,3(2x,F15.8),4X,F5.1)
 80   format (6(1x,e16.10))
 81   format (6(1x,f13.10))
 82   format (8f10.6)
 500  format (' # FRAGMENT ',I2,' FORMED BY ',I3,' BASINS:',
     & 15I4,/,1000(' # ',33x,15I4,/))
 501  format (1x,69('-'))
 502  format (a,'(MinElec,MaxElec)',14x,':',2(1x,I4),
     &        ' (CORE ELECTRONS EXCLUDED)')
 504  format (1x,'#',/,1x,
     &  '# MAXIMUM AND MINIMUM ELECTRON POPULATIONS OF FRAGMENTS')
 503  format (1x,'Cartesian coordinates (bohr) and nuclear charges')
 28   format (//,' # NUMBER OF DETERMINANTS ACTUALLY USED : ',I6)
 29   format (   ' # THERE ARE ',I3,' CORE    ORBITALS :',200(1x,I3))
 31   format (   ' # THERE ARE ',I3,' VALENCE ORBITALS :',200(1x,I3))
 432  format (1x,'#',/,
     &   1x,'# Atom number ',I3,' of group ',I2,1x,
     &  'is equal to atom number ',I3,' of group ',I2,/,1x,'#')
 15   format (/1x,'# Atomic overlap matrix not read for atom ',I3,/)
 16   format (1x,'# Reading AOM for atom ',I3)
 17   format (1x,'# Fatal error reading aom() array:',/,
     &        1x,'# Trying to read aom() for atom ',I3,
     &        ', and there are only ',I3,' atoms')
 21   format (1x,'# Number of electrons               = ',I10,/,
     &        1x,'# Number of molecular orbitals      = ',I10)
 112  format (1x,'#',/,1x,'# ',80('+'),/,
     &        1x,'# ISOPYCNIC LOCALIZATION ROUTINE (GROUPS)',/,1x,'#')
 1121 format (1x,'#',/,1x,'# ',80('+'),/,
     &        1x,'# ISOPYCNIC LOCALIZATION ROUTINE (CENTERS)',/,1x,'#')
 444  format (6(2x,F12.6))
 4444 format (5(1x,F25.18))
 40   format (1x,'#',/,1x,'#',80('+'),/,
     &        1x,'# INVERSE OF THE TRANSPOSE OF THE ROTATION MATRIX',
     &        1x,'FROM NATURAL TO LOCALIZED MOs')
 41   format (//' # ROTATION MATRIX (U) FROM LOCALIZED MOs TO',
     &        ' CANONICAL MOs',/,' # canMO_i =  SUM_j locMO_j U_ji')
 30   format (10(1x,F12.6))
 331  format (' # Canonical MO number ',I3)
 3510 format (1x,"#",/,1x,"# Writing file ","'",a,"' with",
     & " AOMs over Localized MOs")
 3511 format (1x,"#",/,1x,"# Writing file ","'",a,"' with",
     & " AOMs over ",a," Localized MOs")
 351  format (1x,"#",/,1x,"# Writing file ","'",a,"' with",
     & " Localized MOs instead of Natural MOs")
 335  format (1x,'# AOM is renormalized, the reference atom is ',I5)
 221  format (1x,'#',/,1x,'# Automatic analysis of Localized CORE MOs',
     &      /,1x,'# Overlap criterium of Localization = ',F16.10)
 222    format (1x,'# FRAGMENT ',I3,' has ',I3,
     &             ' Localized CORE MOs: ',40I3)
 225    format (1x,'# <LOC_',I3,'|LOC_',I3,'>_',I2,' = ',F16.10)
 1130 format (3(1x,'#',/),1x,'# ',80('+'),/,1x,
     & '# EXACT CALCULATION OF PROBABILITIES WITH CANONICAL MOs',/,
     & ' # (MinRandom,MaxRandom) = ',2(1x,F12.5),/,1x,'#')
 1131 format (3(1x,'#',/),1x,'# ',80('+'),/,1x,
     & '# EXACT CALCULATION OF PROBABILITIES WITH LOCALIZED MOs',/,1x,
     & '# LOCALIZED MOs ARE NOT ORTHOGONAL ==> SUM_i PROB_i .NE. 1',/,
     & ' # (MinRandom,MaxRandom) = ',2(1x,F12.5),/,1x,'#')
 87   format (1x,
     & '# Cutting the WaveFunction expansion in terms of Slater dets',/,
     & 1x,'# Coefficients smaller than ',E17.10,' are neglected',/,
     & 1x,'# Number of Slater determinants before cutting = ',I6,/,
     & 1x,'# Number of Slater determinants after  cutting = ',I6)
 22   format (' # Probability = ',F12.6,4x,'for RSRS ',20I4)
 23   format (' # Approximate EDF assuming that electron populations')
 116  format (1x,'# PROB = ',E17.10,' for the ',a,
     & ' SPINLESS RESONANT STRUCTURE')
 544  format (' # !! AOM is not Ok !! : SUM_n AOM(n,',2I3,') = ',E15.8)
 334  format (/,' # Second & fourth moment orbital spread of ',a)
 10   format (' # MOs are orthogonalized using the OFMO method')
 100  format (' # edf.f: File ',a,' NOT FOUND')
 641  format (
     & " # AOM file with core MOs suppresed is not written",/,
     & " # WFN file with core MOs suppresed is not written",/,
     & " # To write these files, increase the value of MMLINE ",/,
     & " # in the main program 'edf.f' and in the routine 'coremos.f'")
 101  format (' # edf.f: AOMNORM can not be used in UHF calculations',
     & /,' # edf.f: aomnorm ---> false')
 223  format (//1x,'# ',/' # ',I3,' Valence MOs',/,100(1x,'#',20I3,/))
 383  format (/,' # Too many groups (ngroup > ncent',/)
 1133 format (' # USING TOPMOD ',a,' BASINS')
 76   format (///,
     &   ' # This version of EDF only allows "nbasins" and "natoms"',/,
     &   ' # to both be equal to "ncent", the variable read from the',/,
     &   ' # wave function file. In case this is not so, it is',/,
     &   ' # necessary to add ghost atoms in that file. Moreover,',/,
     &   ' # the variables "nbasins" and "natoms", read in from ',/,
     &   ' # the first line of the AOM file must be equal.',//)
 1117 format (' # edf.f: Probabilities are not possible in CCWFNs')
 322  format (1x,
     & '# edf.f: Only OQSEOS order possible in WFNs with only NATMOs')
 985  format (' #',/,
     & ' # +++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' # Effective Oxidation States (EOS) of the fragments',/,
     & ' # +++++++++++++++++++++++++++++++++++++++++++++++++')
 986  format (' # Fragment ',I2,'  EOS = ',F15.8)
 987    format (' # Total            = ',F15.8)
      end
