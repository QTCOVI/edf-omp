c
c-----------------------------------------------------------------------
c
      module space_for_wfxedf
      save

      integer (kind=4)                              :: nrrec
      integer (kind=4)                              :: nedf
      integer (kind=4), allocatable,dimension (:)   :: itypedf
      integer (kind=4), allocatable,dimension (:)   :: icenedf
      real    (kind=8), allocatable,dimension (:)   :: expedf
      real    (kind=8), allocatable,dimension (:)   :: cedf
      end module space_for_wfxedf
c
c-----------------------------------------------------------------------
c
      subroutine allocate_edf_space ()
      USE        space_for_wfxedf
c
      if (.not.allocated(expedf))  allocate (expedf(nedf))
      if (.not.allocated(cedf))    allocate (cedf(nedf))
      if (.not.allocated(itypedf)) allocate (itypedf(nedf))
      if (.not.allocated(icenedf)) allocate (icenedf(nedf))
c
      end subroutine allocate_edf_space
c
c-----------------------------------------------------------------------
c
      subroutine deallocate_edf_space ()
      USE        space_for_wfxedf
c
      if (allocated(expedf))  deallocate (expedf)
      if (allocated(cedf))    deallocate (cedf)
      if (allocated(itypedf)) deallocate (itypedf)
      if (allocated(icenedf)) deallocate (icenedf)
c
      end subroutine deallocate_edf_space
c
c-----------------------------------------------------------------------
c
c-----Read the Electron Density Function (EDF) associated to the core 
c     electrons using the same format employed in the AIMALL code to
c     read the full .wfx file.
c
c-----------------------------------------------------------------------
c
*     subroutine readedf (iwfn,nrrec,thereisEDF)
      subroutine readedf (iwfn,thereisEDF)
      USE space_for_wfxedf
c
      implicit none
      integer (kind=4) mline
      parameter (mline=132)
      integer (kind=4) iwfn,j
      character(len = 50) tag,tag1
      character(len = mline) line,uppcase
      logical thereisEDF
c
c-----------------------------------------------------------------------
c
      tag  = '<ADDITIONAL ELECTRON DENSITY FUNCTION (EDF)>      '
      tag1 = '</ADDITIONAL ELECTRON DENSITY FUNCTION (EDF)>     '
      read (iwfn,'(a)') line
      if (trim(uppcase(line)).eq.trim(tag)) then

        read (iwfn,'(a)') line
        read (iwfn,*,err=1) nedf
        read (iwfn,'(a)') line

        call allocate_edf_space ()

        read (iwfn,'(a)') line
        read (iwfn,*,err=2) (icenedf(j),j=1,nedf)
        read (iwfn,'(a)') line

        read (iwfn,'(a)') line
        read (iwfn,*,err=3) (itypedf(j),j=1,nedf)
        do j=1,nedf
          if (itypedf(j).ne.1) then
            stop '# wfxedf.f: ITYPEDF() # from 1 not allowed yet'
          endif
        enddo
        read (iwfn,'(a)') line
 
        read (iwfn,'(a)') line
        read (iwfn,*,err=4) (expedf(j),j=1,nedf)
        read (iwfn,'(a)') line

        read (iwfn,'(a)') line
        read (iwfn,*,err=5) (cedf(j),j=1,nedf)
        read (iwfn,'(a)') line
        read (iwfn,'(a)') line
        thereisEDF = .true.
c
c-------Determine the number of records of logical unit iwfn associated
c       to the reading of the CORE EDF.
c
        rewind (iwfn)
        do
          read (iwfn,'(a)') line
          if (trim(uppcase(line)).eq.trim(tag)) then
            nrrec=1
          elseif (trim(uppcase(line)).eq.trim(tag1)) then
            nrrec=nrrec+1
            exit
          else
            nrrec=nrrec+1
          endif
        enddo
      else
        thereisEDF = .false.
        backspace (iwfn)
      endif
      return
c
c-----Error in reading data
c
 1    stop ' # wfxedf.f: Bad value read for the NEDF variable'
 2    stop ' # wfxedf.f: Bad values read for ICENEDF()'
 3    stop ' # wfxedf.f: Bad values read for ITYPEDF()'
 4    stop ' # wfxedf.f: Bad values read for EXPEDF()'
 5    stop ' # wfxedf.f: Bad values read for CEDF()'
c
      end
c
c-----------------------------------------------------------------------
c
c-----Print the function defining the core electron density
c
      subroutine writeedf (uout,thereisEDF)
      USE space_for_wfxedf

      implicit none
      real(kind=8), parameter :: pi=3.141592653589793238462643d0
      integer (kind=4) uout,i
      real    (kind=8) corelectrons
      logical thereisEDF,expgtzero
c
      expgtzero = .false.
      do i=1,nedf
        if (itypedf(i) /= 1) expgtzero = .true.
      enddo
      if (thereisEDF) then
        write (uout,900) nedf
        corelectrons = 0.0d+00
        do i=1,nedf
          write (uout,901) i,icenedf(i),itypedf(i),cedf(i),expedf(i)
          corelectrons = corelectrons + cedf(i)*(pi/expedf(i))**1.5d+00
        enddo
        if (.not.expgtzero) write (uout,902) corelectrons
      endif
      return
 900  format (' #',/,
     &  ' # This molecule has a CORE Electron Density Function of ',
     &  I3, ' terms',/,' # Term  Center   Type',13x,'Coef',21x,'Exp')
 901  format (2x,i4,3x,i4,3x,i4,3x,E22.12,3x,E22.12)
 902  format (' # Total number of core electrons = ',12x,E22.12)
      end
