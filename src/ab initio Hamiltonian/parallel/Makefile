# Makefile to compile the subroutines for obtaining the ABINIT wave functions, densities, potentials
# and putting them into form compatible with AI2NBSE codes
#
.SUFFIXES: .f .x .f90
#FC = gfortran

# ---- MPI ----
FC = mpiifort
POPT = -g -traceback -fp-stack-check -check bounds
#
# Executable files from the program files
PREPFILES   = whom0.x datagen.x setup2.x \
	      conugtoux.x orthog.x condens2.x ll.x bridgelad.x bridgegw.x

BSEFILES = spectrum.x zspectrum.x psmnewma_nherm_pro.x \
	
# Object files
#
# --- PREP STAGE ----
DATAGEN_OBJS = datagen.f90 lxb_getval.o getqefiles.o cross.o getrecvec.o 
SETUP_OBJS = setup2.f90 lxb_getval.o enkread.o melread.o cross.o getrecvec.o 
CONVUG_UX_OBJS = conugtoux.f90 lxb_getval.o formx.o gentoreal.o facpowfind.o optim.o chkfftreal.o chkfftrecp.o \
		 jlmfft.o 
ORTHOG_OBJS = orthog.f90 lxb_getval.o dptest.o normalize.o
CONDENS2_OBJS = condens2.f90 lxb_getval.o rho.o solidca.o brcapper.o 
LL_OBJS = ll.f levlou.o lxb_getval.o 
BRIDGELAD_OBJS = bridgelad.f90 lxb_getval.o formx.o cross.o getspacings.o dist.o dublinterp.o
BRIDGEGW_OBJS = bridgegw.f90 lxb_getval.o

# --- BSE STAGE ---
BSE_OBJS = lxb_getval.o formx.o kmapr.o getwgt.o getumat.o bdallck.o uphase.o act.o \
	     act_1e.o smbubble.o smkxc.o smlsauber.o adump.o invdrv.o cinvert.o coulmult.o aclip.o bfnorm.o \
	     getband.o wvfdot.o ctrdfft.o elsdch.o jlmfft.o gradorth.o \
	     inner.o soldrv.o subspace.o velmul.o velmuls.o  
SMNEWMA_OBJS = smnewma.f90 lxb_getval.o formx.o kmapr.o getwgt.o getumat.o bdallck.o uphase.o haydock.o act.o \
	     act_1e.o smbubble.o smkxc.o smlsauber.o adump.o invdrv.o cinvert.o coulmult.o aclip.o bfnorm.o \
	     getband.o wvfdot.o ctrdfft.o elsdch.o jlmfft.o gradorth.o \
	     inner.o soldrv.o subspace.o velmul.o velmuls.o  

# Libraries
LIBS1 = -llapack -lblas
LIBS2 = -lm -lgfortran -lscalapack -llapack -lblas #-llapack 

# Edit this for your installation of fftw3
#FFTI = -I/usr/include/
#FFTWL = -L/usr/lib/x86_64-linux-gnu/ -lfftw3
FFTI = -I/usr/local/almalinux9/inteloneapi/2023.0.0/mkl/2023.0.0/include/fftw/
FFTWL = -lfftw3

# Edit this to link to your stdlib from your gfortran installation
DATAGEN_LIB = #/home/ishiakasem/stdlib/build/gfortran_6C1C75E19A159722/stdlib/libstdlib.a
DATAGEN_INCL = #-I/home/ishiakasem/stdlib/build/gfortran_020AF3609BFD7C6D

# Make and clean
RM = rm -f

all: $(PREPFILES) $(BSEFILES)

clean:
	$(RM) *.o

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

# Compiling executables
#
## ----- PREP STAGE -----

whom0.x: whom0.f
	$(FC) $(OPTIONS) whom0.f -o whom0.x

datagen.x: $(DATAGEN_OBJS)
	$(FC) $(DATAGEN_INCL) -o datagen.x $(DATAGEN_OBJS) $(DATAGEN_LIB) 

setup2.x: $(SETUP_OBJS)
	$(FC) $(OPTIONS) $(SETUP_OBJS) -o setup2.x

conugtoux.x: $(CONVUG_UX_OBJS)
	$(FC) $(OPTIONS) $(CONVUG_UX_OBJS) -o conugtoux.x

orthog.x: $(ORTHOG_OBJS)
	$(FC) $(OPTIONS) $(ORTHOG_OBJS) -o orthog.x

condens2.x: $(CONDENS2_OBJS)
	$(FC) $(OPTIONS) $(CONDENS2_OBJS) -o condens2.x

ll.x: $(LL_OBJS)
	$(FC) $(OPTIONS) $(LL_OBJS) -o ll.x

bridgelad.x: $(BRIDGELAD_OBJS)
	$(FC) $(OPTIONS) $(BRIDGELAD_OBJS) -o bridgelad.x

bridgegw.x: $(BRIDGEGW_OBJS)
	$(FC) $(OPTIONS) $(BRIDGEGW_OBJS) -o bridgegw.x

#---- BSE STAGE -----
smnewma.x: $(SMNEWMA_OBJS)
	$(FC) $(OPTIONS) $(SMNEWMA_OBJS) -o smnewma.x $(LIBS1)

spectrum.x: spectrum.o lxb_getval.o
	$(FC) $(OPTIONS) spectrum.f lxb_getval.o -o spectrum.x

zspectrum.x: zspectrum.o lxb_getval.o
	$(FC) $(OPTIONS) zspectrum.f lxb_getval.o -o zspectrum.x

smnewma_nherm.x: smnewma_nherm.f90 $(BSE_OBJS)
	$(FC) -o smnewma_nherm.x $(OPTIONS) smnewma_nherm.f90 $(BSE_OBJS) $(LIBS1)

# Parallel code with ScaLAPACK

psmnewma_nherm_pro.x: psmnewma_nherm_pro.f90 $(BSE_OBJS)
	$(FC) -o psmnewma_nherm_pro.x $(POPT) psmnewma_nherm_pro.f90 $(BSE_OBJS) $(LIBS2)

# Generating object files
#
#  ooooo PREP STAGE  ooooo
datagen.o: datagen.f90
	$(FC) $(OPTIONS) datagen.f90 -c datagen.o

setup2.o: setup2.f90
	$(FC) $(OPTIONS) setup2.f90 -c setup2.o

conugtoux.o: conugtoux.f90
	$(FC) $(OPTIONS) conugtoux.f90 -c conugtoux.o

orthog.o: orthog.f90
	$(FC) $(OPTIONS) orthog.f90 -c orthog.o

condens2.o: condens2.f90
	$(FC) $(OPTIONS) condens2.f90 -c condens2.o

getqefiles.o: getqefiles.f90
	$(FC) $(OPTIONS) getqefiles.f90 -c getqefiles.o

getrecvec.o: getrecvec.f90
	$(FC) $(OPTIONS) getrecvec.f90 -c getrecvec.o 

cross.o: cross.f90
	$(FC) $(OPTIONS) cross.f90 -c cross.o

lxb_getval.o: lxb_getval.f90
	$(FC) $(OPTIONS) lxb_getval.f90 -c -o lxb_getval.o 

enkread.o: enkread.f90
	$(FC) $(OPTIONS) enkread.f90 -c enkread.o 

melread.o: melread.f90
	$(FC) $(OPTIONS) melread.f90 -c melread.o 

gentoreal.o: gentoreal.f90
	$(FC) $(OPTIONS) gentoreal.f90 -c gentoreal.o 

facpowfind.o: facpowfind.f90
	$(FC) $(OPTIONS) facpowfind.f90 -c facpowfind.o 

optim.o: optim.f90
	$(FC) $(OPTIONS) optim.f90 -c optim.o

chkfftreal.o: chkfftreal.f90
	$(FC) $(OPTIONS) chkfftreal.f90 -c chkfftreal.o 

chkfftrecp.o: chkfftrecp.f90
	$(FC) $(OPTIONS) chkfftrecp.f90 -c chkfftrecp.o 

jlmfft.o: jlmfft.f
	$(FC) $(OPTIONS) jlmfft.f -c -o jlmfft.o 

formx.o: formx.f90
	$(FC) $(OPTIONS) formx.f90 -c -o formx.o 

dptest.o: dptest.f90
	$(FC) $(OPTIONS) dptest.f90 -c dptest.o 

normalize.o: normalize.f90
	$(FC) $(OPTIONS) normalize.f90 -c normalize.o 

rho.o: rho.f90
	$(FC) $(OPTIONS) rho.f90 -c rho.o 

brcapper.o: brcapper.f90
	$(FC) $(OPTIONS) brcapper.f90 -c brcapper.o 

solidca.o: solidca.f
	$(FC) $(OPTIONS) solidca.f -c solidca.o 

levlou.o: levlou.f
	$(FC) $(OPTIONS) levlou.f -c levlou.o 

dublinterp.o: dublinterp.f90
	$(FC) $(OPTIONS) dublinterp.f90 -c dublinterp.o 

dist.o: dist.f90
	$(FC) $(OPTIONS) dist.f90 -c dist.o 

getspacings.o: getspacings.f90
	$(FC) $(OPTIONS) getspacings.f90 -c getspacings.o 

# oooo BSE STAGE ooooo
# need -c -o to run with mpiifort
#
kmapr.o: kmapr.f90
	$(FC) $(OPTIONS) kmapr.f90 -c -o kmapr.o 

getwgt.o: getwgt.f90
	$(FC) $(OPTIONS) getwgt.f90 -c -o getwgt.o 

getumat.o: getumat.f90
	$(FC) $(OPTIONS) getumat.f90 -c -o getumat.o 

bdallck.o: bdallck.f90
	$(FC) $(OPTIONS) bdallck.f90 -c -o bdallck.o 

uphase.o: uphase.f90
	$(FC) $(OPTIONS) uphase.f90 -c -o uphase.o 

haydock.o: haydock.f90
	$(FC) $(OPTIONS) haydock.f90 -c -o haydock.o 

act.o: act.f90
	$(FC) $(OPTIONS) act.f90 -c -o act.o 

act_1e.o: act_1e.f90
	$(FC) $(OPTIONS) act_1e.f90 -c -o act_1e.o 

smbubble.o: smbubble.f90
	$(FC) $(OPTIONS) smbubble.f90 -c -o smbubble.o 

smkxc.o: smkxc.f90
	$(FC) $(OPTIONS) smkxc.f90 -c -o smkxc.o 

smlsauber.o: smlsauber.f90
	$(FC) $(OPTIONS) smlsauber.f90 -c -o smlsauber.o 

adump.o: adump.f90
	$(FC) $(OPTIONS) adump.f90 -c -o adump.o 

invdrv.o: invdrv.f90
	$(FC) $(OPTIONS) invdrv.f90 -c -o invdrv.o 

cinvert.o: cinvert.f
	$(FC) $(OPTIONS) cinvert.f -c -o cinvert.o 

coulmult.o: coulmult.f90
	$(FC) $(OPTIONS) coulmult.f90 -c -o coulmult.o 

aclip.o: aclip.f90
	$(FC) $(OPTIONS) aclip.f90 -c -o aclip.o 

bfnorm.o: bfnorm.f90
	$(FC) $(OPTIONS) bfnorm.f90 -c -o bfnorm.o 

getband.o: getband.f90
	$(FC) $(OPTIONS) getband.f90 -c -o getband.o 

ldriver.o: ldriver.f90
	$(FC) $(OPTIONS) ldriver.f90 -c -o ldriver.o 

lpart1.o: lpart1.f90
	$(FC) $(OPTIONS) lpart1.f90 -c -o lpart1.o 

lpart2.o: lpart2.f90
	$(FC) $(OPTIONS) lpart2.f90 -c -o lpart2.o 

lpart3.o: lpart3.f90
	$(FC) $(OPTIONS) lpart3.f90 -c -o lpart3.o 

tid.o: tid.f90
	$(FC) $(OPTIONS) tid.f90 -c -o tid.o 

wvfdot.o: wvfdot.f90
	$(FC) $(OPTIONS) wvfdot.f90 -c -o wvfdot.o 

ctrdfft.o: ctrdfft.f90
	$(FC) $(OPTIONS) ctrdfft.f90 -c -o ctrdfft.o 

elsdch.o: elsdch.f
	$(FC) $(OPTIONS) elsdch.f -c -o elsdch.o 

gradorth.o: gradorth.f
	$(FC) $(OPTIONS) gradorth.f -c -o gradorth.o 

inner.o: inner.f
	$(FC) $(OPTIONS) inner.f -c -o inner.o 

soldrv.o: soldrv.f
	$(FC) $(OPTIONS) soldrv.f -c -o soldrv.o 

subspace.o: subspace.f
	$(FC) $(OPTIONS) subspace.f -c -o subspace.o 

velmul.o: velmul.f90
	$(FC) $(OPTIONS) velmul.f90 -c -o velmul.o 

velmuls.o: velmuls.f90
	$(FC) $(OPTIONS) velmuls.f90 -c -o velmuls.o 

spectrum.o: spectrum.f
	$(FC) $(OPTIONS) spectrum.f -c -o spectrum.o 

zspectrum.o: zspectrum.f
	$(FC) $(OPTIONS) zspectrum.f -c -o zspectrum.o 
clean:
	rm -f *.x *.o
