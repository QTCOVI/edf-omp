#! /usr/bin/make

# ++++++++++++++++++++++++++++   gfortran compilation ++++++++++++++++++++++++++++++++++++++
# 
# Static compilation (Then, openmp CAN NOT be used)
#
   obj    = edf
   all : $(obj)
#  FCOMPL = gfortran -fopenmp -Wall -Wunused-label -Wmaybe-uninitialized -fbounds-check
   FCOMPL = gfortran -fopenmp -fbounds-check
   CCOMPL = gcc 
   CFLAGC =
   opt    = -O
   SYSDEP =
   ldflag = -static -static-libgcc /usr/lib/x86_64-linux-gnu/lapack/liblapack.a /usr/lib/x86_64-linux-gnu/blas/libblas.a
#
#------------------------------------------------------------------------------------------------
# 
# Non-Static compilation (Then, openmp CAN be used)
#
# obj    = edf
# all : $(obj)
# FCOMPL = gfortran -fopenmp
# CCOMPL = gcc -fopenmp
# CFLAGC =
# opt    = -O
# SYSDEP =
# ldflag = /usr/lib/x86_64-linux-gnu/lapack/liblapack.a /usr/lib/x86_64-linux-gnu/blas/libblas.a
# 
# ++++++++++++++++++++++++++++   gfortran compilation ++++++++++++++++++++++++++++++++++++++







# Implicit rules

.SUFFIXES:
.SUFFIXES: .vec .f .c .s .o .f90

.f.o:
	$(FCOMPL) -c $(FFLAGC) $(opt) -o $@ $<

.f90.o:
	$(FCOMPL) -c $(FFLAGC) $(opt) -o $@ $<

.s.o:
	as -o $@ $<

.c.o:
	$(CCOMPL) -c $(CFLAGC) -o $@ $<


# Object files

OBJECTS = mod_sym.o mod_periodic.o symarea.o wfnbasis.o wfncoef.o rdm1space.o datatm.o \
          primgto.o cidetarea.o spaceconf.o rsspace.o bondarea.o dafhspace.o sgarea.o wfxedf.o \
          aomspin.o cdafhmos.o wrtwfn.o qqsort.o qcksort.o lmowfn.o jacobi.o \
          cutwfn.o srhon.o colord.o edf.o edfstd.o edfcrit.o edfx.o error.o leng.o nbasab.o \
          xnbasab.o rnprobs.o random.o rdwfn.o semilla.o atoi.o isdigit.o \
          setdble.o setint.o setword.o sndelta.o locdeloc.o ndelta.o binedf.o  \
          calcedfd.o xcalcedfd.o rcalcedf.o xcalcedf.o timer.o \
          uppcase.o iqcksort.o ruedmis.o getdate.o aomlowdin.o aomulliken.o \
          aomindef.o sabxyz.o locsome.o suppcas.o suppress.o delmo.o isopycore.o  nlmfill.o \
          sdwrsrs.o wrsrs.o probres.o configs.o  sdwsumdet.o sumdetp.o \
          twocendi.o cindgr.o ctwocen.o lengrec.o \
          compedfstd.o xcritcplx.o critcplx.o cmplxedf.o readcmplx.o readaomstd.o readaom.o znetlib.o gendafhstd.o gendafh.o \
          dafhstd.o twobascrit.o dafhcrit.o dafh.o dnostd.o dnocrit.o xdnoedf.o dnoedf.o cmplxlo.o locastd.o critlocmplx.o locmplx.o cgedfp.o rhocond.o allrhoc.o rhoalls.o \
          indepgr.o nbasabw.o optedf.o  praxis.o corredf.o hfedf.o gedfit.o \
          wfnmax.o fdlbfgs.o lbfgs.o molorb0.o molorb1.o molorb2.o mxwfn.o \
          gtogto.o etijcalc.o chsort.o searchstring.o computerdm.o \
          aomtest.o canlmo.o aomindefrho.o aombecke.o sphere_lebedev_rule.o calcmos.o \
          multipole.o newdrvmult.o uhfchemloc.o chemloc.o itdafhc.o sdwlewis.o ulewisoqs.o lewisoqs.o casoqsloc.o oqsloc.o itdafh.o \
          uhfitaomloc.o itaomloc.o \
          confrag.o connect.o promrhow.o netrhow.o mindefw.o atomicrho.o \
          normalizemos.o ofmortho.o errdsyev.o eiganalysis.o lower.o  newsym.o \
          sym.o pntgrp.o symgr.o typeop.o eispack.o transor.o orbfun.o auxfile.o \
          dosym.o compdafh.o fourchar.o coremos.o mutent.o timestamp.o \
          sopenrho.o openrho.o dafhdrv.o dafhdo.o topoint.o drvother.o udrvother.o otheroqs.o \
          detlapack.o lulapack.o ngen.o casuhfoqsloc.o uhfoqsloc.o casuhfdrvloc.o  uhfdrvloc.o \
          uhfnorm.o fnofrag.o neweos.o sdweos.o nfnoeos.o fnoeos.o lubksb.o ludcmp.o mylapack.o variance.o \
          diagcov.o sdwdet.o calcres.o gengen.o aomnorma.o dirprob.o pfromab.o \
          locrohf.o writeaom.o testaom.o wbec.o wscuseria.o

# Clean objects

clean:
	@rm -f core $(VECTOR) $(OBJECTS) *.mod $(obj) edf-omp

# Main program target

$(obj): $(SYSDEP) $(VECTOR) $(OBJECTS)
	$(FCOMPL) $(FFLAGC) -o $(obj) $(OBJECTS) $(SYSDEP) $(ldflag)

# OpenMP-enabled dynamic-link target.
# Drops -static/-static-libgcc so libgomp resolves dynamically and the
# user-level !$omp directives actually run multi-threaded.
# Runtime: set OMP_NUM_THREADS=<cores>. If linking later against a threaded
# BLAS (OpenBLAS/MKL) set OPENBLAS_NUM_THREADS=1 (or equivalent) to avoid
# over-subscription when hot loops are already threaded at the Fortran level.
edf-omp: $(SYSDEP) $(VECTOR) $(OBJECTS)
	$(FCOMPL) $(FFLAGC) -fopenmp -o edf-omp $(OBJECTS) $(SYSDEP) \
	    /usr/lib/x86_64-linux-gnu/lapack/liblapack.a \
	    /usr/lib/x86_64-linux-gnu/blas/libblas.a

# Include files dependencies

          allrhoc.o:   implicit.inc
          allrhoc.o:        wfn.inc
          allrhoc.o:       corr.inc

         aombecke.o:        wfn.inc
         aombecke.o:     datatm.inc
         aomindef.o:   implicit.inc
         aomindef.o:        wfn.inc
         aomindef.o:      param.inc
         aomindef.o:  constants.inc

      aomindefrho.o:        wfn.inc
      aomindefrho.o:     datatm.inc

        aomlowdin.o:   implicit.inc
        aomlowdin.o:        wfn.inc
        aomlowdin.o:      param.inc
        aomlowdin.o:  constants.inc

          aomtest.o:   implicit.inc
          aomtest.o:  constants.inc
          aomtest.o:      error.inc

       aomulliken.o:   implicit.inc
       aomulliken.o:        wfn.inc
       aomulliken.o:      param.inc
       aomulliken.o:  constants.inc

           binedf.o:   implicit.inc
           binedf.o:      param.inc
           binedf.o:  constants.inc
           binedf.o:    lengrec.inc

         calcedfd.o:   implicit.inc
         calcedfd.o:      param.inc
         calcedfd.o:  constants.inc
         calcedfd.o:    lengrec.inc

          calcmos.o:   implicit.inc
          calcmos.o:      param.inc
          calcmos.o:        wfn.inc
          calcmos.o:    primgto.inc

           canlmo.o:   implicit.inc
           canlmo.o:  constants.inc
           canlmo.o:      error.inc

         cdafhmos.o:   implicit.inc
         cdafhmos.o:      param.inc
         cdafhmos.o:        wfn.inc
         cdafhmos.o:       corr.inc

           cgedfp.o:   implicit.inc
           cgedfp.o:      param.inc
           cgedfp.o:  constants.inc
           cgedfp.o:        wfn.inc
           cgedfp.o:       corr.inc

          chemloc.o:   implicit.inc
          chemloc.o:  constants.inc
          chemloc.o:      error.inc

          uhfchemloc.o:   implicit.inc
          uhfchemloc.o:  constants.inc
          uhfchemloc.o:      error.inc
          uhfchemloc.o:      wfn.inc
          uhfchemloc.o:      corr.inc

           chsort.o:   implicit.inc
           chsort.o:  constants.inc
           chsort.o:     stderr.inc
           cindgr.o:   implicit.inc
           cindgr.o:      param.inc
           cindgr.o:  constants.inc
           cindgr.o:    lengrec.inc

         cmplxedf.o:   implicit.inc
         cmplxedf.o:      param.inc
         cmplxedf.o:  constants.inc

         critcplx.o:   implicit.inc
         critcplx.o:      param.inc
         critcplx.o:  constants.inc

         xcritcplx.o:   implicit.inc
         xcritcplx.o:      param.inc
         xcritcplx.o:  constants.inc


         compedfstd.o:   implicit.inc
         compedfstd.o:      param.inc
         compedfstd.o:  constants.inc

          cmplxlo.o:   implicit.inc
          cmplxlo.o:  constants.inc
          cmplxlo.o:      error.inc

          configs.o:   implicit.inc
          configs.o:  constants.inc

          coremos.o:   implicit.inc
          coremos.o:      mline.inc
          coremos.o:    lengrec.inc

          corredf.o:   implicit.inc
          corredf.o:      param.inc
          corredf.o:  constants.inc

          ctwocen.o:   implicit.inc
          ctwocen.o:      param.inc
          ctwocen.o:  constants.inc

           cutwfn.o:   implicit.inc
           cutwfn.o:  constants.inc
           cutwfn.o:      param.inc

           dafhdo.o:   implicit.inc
           dafhdo.o:  constants.inc
           dafhdo.o:      param.inc
           dafhdo.o:        wfn.inc
           dafhdo.o:      error.inc

         otheroqs.o:   implicit.inc
         otheroqs.o:  constants.inc
         otheroqs.o:      param.inc
         otheroqs.o:        wfn.inc
         otheroqs.o:      error.inc


             dafh.o:   implicit.inc
             dafh.o:      param.inc
             dafh.o:  constants.inc
             dafh.o:     stderr.inc

             dafhcrit.o:   implicit.inc
             dafhcrit.o:      param.inc
             dafhcrit.o:  constants.inc
             dafhcrit.o:     stderr.inc

             twobascrit.o:   implicit.inc
             twobascrit.o:      param.inc
             twobascrit.o:  constants.inc
             twobascrit.o:     stderr.inc

          dafhstd.o:   implicit.inc
          dafhstd.o:      param.inc
          dafhstd.o:  constants.inc

           datatm.o:   implicit.inc
           datatm.o:     datatm.inc

            delmo.o:   implicit.inc
            delmo.o:      param.inc
            delmo.o:  constants.inc
            delmo.o:        wfn.inc
            delmo.o:       corr.inc
            delmo.o:     stderr.inc

           dnoedf.o:   implicit.inc
           dnoedf.o:      param.inc
           dnoedf.o:  constants.inc
           dnoedf.o:     stderr.inc

           xdnoedf.o:   implicit.inc
           xdnoedf.o:      param.inc
           xdnoedf.o:  constants.inc
           xdnoedf.o:     stderr.inc

           dnocrit.o:   implicit.inc
           dnocrit.o:      param.inc
           dnocrit.o:  constants.inc
           dnocrit.o:     stderr.inc

           dnostd.o:   implicit.inc
           dnostd.o:      param.inc
           dnostd.o:  constants.inc

          dafhdrv.o:   implicit.inc
          dafhdrv.o:        wfn.inc

            dosym.o:     datatm.inc

          drvmult.o:      param.inc

              edf.o:   implicit.inc
              edf.o:      param.inc
              edf.o:        wfn.inc
              edf.o:       corr.inc
              edf.o:     stderr.inc
              edf.o:  constants.inc
              edf.o:    lengrec.inc
              edf.o:       fact.inc
              edf.o:        sym.inc
              edf.o:     datatm.inc
              edf.o:      integ.inc

             edfx.o:   implicit.inc
             edfx.o:        wfn.inc
             edfx.o:     stderr.inc
             edfx.o:  constants.inc

             edfcrit.o:   implicit.inc
             edfcrit.o:        wfn.inc
             edfcrit.o:     stderr.inc
             edfcrit.o:  constants.inc

           edfstd.o:   implicit.inc
           edfstd.o:        wfn.inc
           edfstd.o:     stderr.inc
           edfstd.o:  constants.inc

            error.o:     stderr.inc
            error.o:      error.inc

          fdlbfgs.o:   implicit.inc
          fdlbfgs.o:      error.inc

           gedfit.o:   implicit.inc
           gedfit.o:       corr.inc
           gedfit.o:      param.inc
           gedfit.o:  constants.inc

          gendafh.o:   implicit.inc
          gendafh.o:      param.inc
          gendafh.o:  constants.inc
          gendafh.o:     stderr.inc

          gendafhstd.o:   implicit.inc
          gendafhstd.o:      param.inc
          gendafhstd.o:  constants.inc
          gendafhstd.o:     stderr.inc

           gtogto.o:   implicit.inc
           gtogto.o:      param.inc
           gtogto.o:        wfn.inc
           gtogto.o:    primgto.inc

            hfedf.o:   implicit.inc
            hfedf.o:      param.inc
            hfedf.o:  constants.inc

          indepgr.o:   implicit.inc
          indepgr.o:      param.inc

        isopycore.o:   implicit.inc
        isopycore.o:      param.inc
        isopycore.o:        wfn.inc
        isopycore.o:       corr.inc
        isopycore.o:     stderr.inc
        isopycore.o:  constants.inc

         itaomloc.o:   implicit.inc
         itaomloc.o:  constants.inc
         itaomloc.o:      error.inc

         uhfitaomloc.o:   implicit.inc
         uhfitaomloc.o:  constants.inc
         uhfitaomloc.o:      error.inc
         uhfitaomloc.o:      corr.inc
         uhfitaomloc.o:      wfn.inc

          itdafhc.o:   implicit.inc
          itdafhc.o:  constants.inc
          itdafhc.o:      param.inc
          itdafhc.o:        wfn.inc
          itdafhc.o:      error.inc

           oqsloc.o:   implicit.inc
           oqsloc.o:  constants.inc
           oqsloc.o:      param.inc
           oqsloc.o:        wfn.inc
           oqsloc.o:      error.inc

        casoqsloc.o:   implicit.inc
        casoqsloc.o:  constants.inc
        casoqsloc.o:      param.inc
        casoqsloc.o:        wfn.inc
        casoqsloc.o:      error.inc

        lewisoqs.o:   implicit.inc
        lewisoqs.o:  constants.inc
        lewisoqs.o:      param.inc
        lewisoqs.o:        wfn.inc
        lewisoqs.o:       corr.inc
        lewisoqs.o:      error.inc

        ulewisoqs.o:   implicit.inc
        ulewisoqs.o:  constants.inc
        ulewisoqs.o:      param.inc
        ulewisoqs.o:        wfn.inc
        ulewisoqs.o:       corr.inc
        ulewisoqs.o:      error.inc

        sdwlewis.o:   implicit.inc
        sdwlewis.o:  constants.inc
        sdwlewis.o:      param.inc
        sdwlewis.o:        wfn.inc
        sdwlewis.o:       corr.inc
        sdwlewis.o:      error.inc


        uhfoqsloc.o:   implicit.inc
        uhfoqsloc.o:  constants.inc
        uhfoqsloc.o:      param.inc
        uhfoqsloc.o:      error.inc
        uhfoqsloc.o:       corr.inc
        uhfoqsloc.o:      mline.inc
        uhfoqsloc.o:      wfn.inc

        casuhfoqsloc.o:   implicit.inc
        casuhfoqsloc.o:  constants.inc
        casuhfoqsloc.o:      param.inc
        casuhfoqsloc.o:      error.inc
        casuhfoqsloc.o:       corr.inc
        casuhfoqsloc.o:      mline.inc
        casuhfoqsloc.o:      wfn.inc

        uhfdrvloc.o:   implicit.inc
        uhfdrvloc.o:  constants.inc
        uhfdrvloc.o:      param.inc
        uhfdrvloc.o:      error.inc
        uhfdrvloc.o:       corr.inc
        uhfdrvloc.o:      mline.inc

        casuhfdrvloc.o:   implicit.inc
        casuhfdrvloc.o:  constants.inc
        casuhfdrvloc.o:      param.inc
        casuhfdrvloc.o:      error.inc
        casuhfdrvloc.o:       corr.inc
        casuhfdrvloc.o:      mline.inc

        udrvother.o:   implicit.inc
        udrvother.o:  constants.inc
        udrvother.o:      param.inc
        udrvother.o:      error.inc
        udrvother.o:       corr.inc
        udrvother.o:      mline.inc

        drvother.o:   implicit.inc
        drvother.o:  constants.inc
        drvother.o:      param.inc
        drvother.o:      error.inc
        drvother.o:       corr.inc
        drvother.o:      mline.inc

           itdafh.o:   implicit.inc
           itdafh.o:  constants.inc
           itdafh.o:      param.inc
           itdafh.o:        wfn.inc
           itdafh.o:      error.inc

          lengrec.o:    lengrec.inc

           lmowfn.o:   implicit.inc
           lmowfn.o:      param.inc
           lmowfn.o:  constants.inc
           lmowfn.o:    lengrec.inc

          locmplx.o:   implicit.inc
          locmplx.o:  constants.inc
          locmplx.o:      error.inc
          locmplx.o:     stderr.inc

          critlocmplx.o:   implicit.inc
          critlocmplx.o:  constants.inc
          critlocmplx.o:      error.inc
          critlocmplx.o:     stderr.inc

          locastd.o:   implicit.inc
          locastd.o:  constants.inc
          locastd.o:      error.inc
          locastd.o:     stderr.inc

          locsome.o:   implicit.inc
          locsome.o:      param.inc
          locsome.o:  constants.inc
          locsome.o:        wfn.inc
          locsome.o:       corr.inc
          locsome.o:     stderr.inc

          mindefw.o:   implicit.inc
          mindefw.o:      param.inc
          mindefw.o:        wfn.inc
          mindefw.o:    primgto.inc

          molorb0.o:   implicit.inc
          molorb0.o:      param.inc
          molorb0.o:  constants.inc
          molorb0.o:        wfn.inc

          molorb1.o:   implicit.inc
          molorb1.o:      param.inc
          molorb1.o:  constants.inc
          molorb1.o:        wfn.inc

          molorb2.o:   implicit.inc
          molorb2.o:      param.inc
          molorb2.o:  constants.inc
          molorb2.o:        wfn.inc

        multipole.o:   implicit.inc
        multipole.o:    primgto.inc

           mutent.o:   implicit.inc

            mxwfn.o:   implicit.inc
            mxwfn.o:      param.inc
            mxwfn.o:  constants.inc
            mxwfn.o:        opt.inc
            mxwfn.o:       corr.inc
            mxwfn.o:        wfn.inc

           nbasab.o:   implicit.inc
           nbasab.o:      param.inc
           nbasab.o:  constants.inc

          nbasabw.o:   implicit.inc
          nbasabw.o:      param.inc
          nbasabw.o:  constants.inc

           ndelta.o:   implicit.inc
           ndelta.o:      param.inc
           ndelta.o:  constants.inc
           ndelta.o:   implicit.inc
           ndelta.o:      param.inc
           ndelta.o:  constants.inc

          netrhow.o:   implicit.inc
          netrhow.o:      param.inc
          netrhow.o:        wfn.inc
          netrhow.o:    primgto.inc

       newdrvmult.o:      param.inc

          nlmfill.o:      param.inc

          openrho.o:   implicit.inc
          openrho.o:  constants.inc

          sopenrho.o:   implicit.inc
          sopenrho.o:  constants.inc
          sopenrho.o:      param.inc
          sopenrho.o:       corr.inc
          sopenrho.o:      mline.inc

           optedf.o:   implicit.inc
           optedf.o:      param.inc
           optedf.o:  constants.inc

           orbfun.o:   implicit.inc
           orbfun.o:      param.inc
           orbfun.o:        wfn.inc

           pntgrp.o:   implicit.inc
           pntgrp.o:   implicit.inc
           pntgrp.o:   implicit.inc
           pntgrp.o:   implicit.inc
           pntgrp.o:   implicit.inc
           pntgrp.o:   implicit.inc
           pntgrp.o:   implicit.inc
           pntgrp.o:   implicit.inc
           pntgrp.o:   implicit.inc
           pntgrp.o:   implicit.inc

         promrhow.o:   implicit.inc
         promrhow.o:      param.inc
         promrhow.o:        wfn.inc
         promrhow.o:    primgto.inc

            probres.o:   implicit.inc
            probres.o:       fact.inc
            probres.o:  constants.inc
            probres.o:       fact.inc
            probres.o:  constants.inc

            wrsrs.o:   implicit.inc
            wrsrs.o:       fact.inc
            wrsrs.o:  constants.inc

            gengen.o:  implicit.inc

         sdwrsrs.o:   implicit.inc
         sdwrsrs.o:       fact.inc
         sdwrsrs.o:  constants.inc

          qcksort.o:     stderr.inc

           qqsort.o:   implicit.inc
           qqsort.o:  constants.inc
           qqsort.o:     stderr.inc

           random.o:   implicit.inc
           random.o:  constants.inc

         rcalcedf.o:   implicit.inc
         rcalcedf.o:      param.inc
         rcalcedf.o:  constants.inc
         rcalcedf.o:    lengrec.inc

            rdwfn.o:      param.inc
            rdwfn.o:        wfn.inc
            rdwfn.o:       corr.inc
            rdwfn.o:    lengrec.inc
            rdwfn.o:       fact.inc
            rdwfn.o:    primgto.inc

          readaom.o:  constants.inc

          readaomstd.o:  constants.inc

        readcmplx.o:  constants.inc

          rhoalls.o:   implicit.inc
          rhoalls.o:        wfn.inc
          rhoalls.o:       corr.inc

          rhocond.o:   implicit.inc
          rhocond.o:      param.inc
          rhocond.o:  constants.inc
          rhocond.o:        wfn.inc
          rhocond.o:       corr.inc

          rnprobs.o:   implicit.inc
          rnprobs.o:   implicit.inc

          ruedmis.o:   implicit.inc
          ruedmis.o:  constants.inc
          ruedmis.o:      error.inc

           setint.o:   implicit.inc

          sndelta.o:   implicit.inc
          sndelta.o:      param.inc
          sndelta.o:  constants.inc
          sndelta.o:   implicit.inc
          sndelta.o:      param.inc
          sndelta.o:  constants.inc

            srhon.o:   implicit.inc
            srhon.o:  constants.inc
            srhon.o:      param.inc
            srhon.o:    lengrec.inc

          sdwsumdet.o:   implicit.inc
          sdwsumdet.o:  constants.inc

          sumdetpsdw.o:   implicit.inc
          sumdetpsdw.o:  constants.inc

         suppress.o:   implicit.inc
         suppress.o:      param.inc
         suppress.o:  constants.inc
         suppress.o:        wfn.inc
         suppress.o:       corr.inc
         suppress.o:     stderr.inc

         suppcas.o:   implicit.inc
         suppcas.o:      param.inc
         suppcas.o:  constants.inc
         suppcas.o:        wfn.inc
         suppcas.o:       corr.inc
         suppcas.o:     stderr.inc
         suppcas.o:     lengrec.inc

              sym.o:   implicit.inc
              sym.o:      param.inc
              sym.o:      integ.inc
              sym.o:        wfn.inc
              sym.o:        sym.inc
              sym.o:     datatm.inc
              sym.o:   implicit.inc

            symgr.o:   implicit.inc

            timer.o:      error.inc

          transor.o:   implicit.inc
          transor.o:      param.inc
          transor.o:      integ.inc
          transor.o:        wfn.inc
          transor.o:        sym.inc

         twocendi.o:   implicit.inc
         twocendi.o:      param.inc
         twocendi.o:  constants.inc

           typeop.o:   implicit.inc
           typeop.o:      param.inc
           typeop.o:   implicit.inc
           typeop.o:   implicit.inc
           typeop.o:   implicit.inc
           typeop.o:   implicit.inc
           typeop.o:   implicit.inc
           typeop.o:      param.inc
           typeop.o:   implicit.inc
           typeop.o:      param.inc
           typeop.o:   implicit.inc

           wfnmax.o:   implicit.inc
           wfnmax.o:      param.inc
           wfnmax.o:  constants.inc
           wfnmax.o:        opt.inc
           wfnmax.o:       corr.inc

           wrtwfn.o:   implicit.inc
           wrtwfn.o:      param.inc
           wrtwfn.o:        wfn.inc
           wrtwfn.o:       corr.inc

        xcalcedfd.o:   implicit.inc
        xcalcedfd.o:      param.inc
        xcalcedfd.o:  constants.inc
        xcalcedfd.o:    lengrec.inc

         xcalcedf.o:   implicit.inc
         xcalcedf.o:      param.inc
         xcalcedf.o:  constants.inc
         xcalcedf.o:    lengrec.inc

          xnbasab.o:   implicit.inc
          xnbasab.o:      param.inc
          xnbasab.o:  constants.inc


         aombecke.o:      mline.inc   
         aomindef.o:      mline.inc   
      aomindefrho.o:      mline.inc   
        aomlowdin.o:      mline.inc   
       aomulliken.o:      mline.inc   
        oldbinedf.o:      mline.inc   
           binedf.o:      mline.inc   
         calcedfd.o:      mline.inc   
         cdafhmos.o:      mline.inc
           cgedfp.o:      mline.inc   
           cindgr.o:      mline.inc   
         compdafh.o:      mline.inc   
          configs.o:      mline.inc   
          corredf.o:      mline.inc   
          ctwocen.o:      mline.inc   
           dafhdo.o:      mline.inc   
          dafhdrv.o:      mline.inc   
              edf.o:      mline.inc   
             edfx.o:      mline.inc   
           edfstd.o:      mline.inc   
          indepgr.o:      mline.inc   
          itdafhc.o:      mline.inc   
           itdafh.o:      mline.inc   
        itopenloc.o:      mline.inc   
         otheroqs.o:      mline.inc   
           lmowfn.o:      mline.inc   
          openrho.o:      mline.inc
         sopenrho.o:      mline.inc   
           optedf.o:      mline.inc   
         rcalcedf.o:      mline.inc   
            rdwfn.o:      mline.inc   
          readaom.o:      mline.inc   
       readaomstd.o:      mline.inc   
          rhocond.o:      mline.inc   
         suppress.o:      mline.inc   
         twocendi.o:      mline.inc   
           wfnmax.o:      mline.inc   
           wrtwfn.o:      mline.inc   
        xcalcedfd.o:      mline.inc   
         xcalcedf.o:      mline.inc   

              fnoeos.o:   implicit.inc
              fnoeos.o:   mline.inc

              nfnoeos.o:   implicit.inc
              nfnoeos.o:   mline.inc

              sdweos.o:   implicit.inc
              sdweos.o:   mline.inc

              neweos.o:   implicit.inc
              neweos.o:   mline.inc

              fnofrag.o:   implicit.inc
              fnofrag.o:   mline.inc

              variance.o:   param.inc
              
               diagcov.o:   param.inc

                sdwdet.o: implicit.inc
                sdwdet.o: constants.inc

                aomnorma.o: implicit.inc
                aomnorma.o: wfn.inc
                aomnorma.o: param.inc
                aomnorma.o: corr.inc
 
                wbec.o : param.inc
                wbec.o : wfn.inc
                wbec.o : datatm.inc
 
                wscuseria.o : param.inc
                wscuseria.o : wfn.inc
                wscuseria.o : datatm.inc
