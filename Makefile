#### Makefile for channel code #### 
EXEC = chnl.x
TGZFILE = all.tgz
#
SRC = \
advance.F io.F pstep.F timers.F viscxyz.F courant.F  rhs.F  \
main.F initial.F divg.F partial.F  lagrangian_pts_gen.F \
stats.F tt_rhs.F post_proc.F nltrms_up.F gen_helmholz.F     \
enlred.F fft2d_new.F idminmax.F  \
ibm.F ibm_tt_new.F buildup.F\
#sphere_nonuniform.F del_fn.F ibm.F modified_windows_fn.F \
#
INC = common.inc fft.inc flags.inc global.inc timers.inc 
#
XTRA = Makefile chnl.anz_ini \
       README README_eqns README_gbal README_time README_old
#
OBJS = $(SRC:.F=.o)
#
# Note: Requires "module load intel fftw" prior to "make UFHPC=1 MP=1"
# 
ifdef UFHPC
MKLDIR     = $(HPC_MKL_DIR)
MKLROOT    = $(HPC_MKL_DIR)
MKLLIBDIR  = ${MKLDIR}/lib/intel64
MKLINCDIR  = ${MKLDIR}/include
MKLFFTWINC = ${MKLINCDIR}/fftw
LIBS       = -L${MKLLIBDIR} \
             -Wl,--start-group \
   	     -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
             -Wl,--end-group 
FLAGS      = -mcmodel=medium -shared-intel -fpp -DIFC -DFFTW3 -I${MKLFFTWINC}
#DEBUGFLAGS = -g -C -traceback
DEBUGFLAGS = -traceback
OPTFLAGS   = -O2 -msse3 -axcore-avx2,core-avx-i 
MPFLAGS    = -qopenmp -D OPENMP
FCMP       = ifort
FCSP       = ifort
endif
#
TARFILE    = chnl.tar
#
ifdef MP
  FC = $(FCMP)
  FLAGS += $(MPFLAGS)
else
  FC = $(FCSP)
endif
ifdef DEBUG
  FLAGS += $(DEBUGFLAGS)
else
  FLAGS += $(OPTFLAGS)
endif

%.o: %.F $(INC)
	$(FC) -c $(FLAGS) $<

$(EXEC): $(OBJS)
	$(FC) -o $@ $(FLAGS) $^ $(LIBS)

dist:
	-rm -f $(TARFILE).gz
	tar -cf $(TARFILE) $(SRC) $(INC) $(XTRA)
	gzip $(TARFILE)

get:
	scp mcantero@cee-zze6.cee.uiuc.edu:~/research/spectral_code/front_version/code/code_mixed9.0_aniso/$(TGZFILE) .
	tar -zxvf $(TGZFILE)
	rm -f $(TGZFILE)

clean:
	-rm -f $(OBJS)
