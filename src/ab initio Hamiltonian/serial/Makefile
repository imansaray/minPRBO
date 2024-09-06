# Makefile to compile the subroutines for obtaining the ABINIT wave functions, densities, potentials
# and putting them into form compatible with Eric SHirley's NBSE codes
# When using Erics latest smlsauber, its correct version is in ldriver.f90.
.SUFFIXES: .f .x .f90
FC = gfortran

# Executable files from the program files
PREPFILES   = tmelnew.x frmtwf.x whom0.x nelectron.x rhoofg.x datagen.x setup2.x \
	      conugtoux.x orthog.x condens2.x ll.x bridgelad.x bridgegw.x

BSEFILES = smnewma.x smnewma_nherm.x spectrum.x zspectrum.x
#BSEFILES = smnewma_nherm.x spectrum.x zspectrum.x
#BSEFILES = smnewma.x spectrum.x zspectrum.x
#BSEFILES = smnewma_herm_no.x smnewma_herm_resolv.x spectrum.x zspectrum.x
#BSEFILES = smnewma_herm_resolv.x spectrum.x zspectrum.x

# Scripts, will leave out for now
#SCRIPTS = dendip.script int2.pl bvecs.pl dos.pl gvecs2.pl fermi.pl ofermi.pl dendip.pl

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
BSE_OBJS = lxb_getval.o formx.o kmapr.o getwgt.o getumat.o bdallck.o uphase.o sizereport.o haydock_herm_resolv.o \
	     act.o act_1e.o smbubble.o smkxc.o smlsauber.o adump.o invdrv.o cinvert.o coulmult.o aclip.o bfnorm.o \
	     getband.o wvfdot.o ctrdfft.o elsdch.o jlmfft.o gradorth.o \
	     inner.o soldrv.o subspace.o velmul.o velmuls.o  

#BSE_OBJS = lxb_getval.o formx.o kmapr.o getwgt.o getumat.o bdallck.o uphase.o sizereport.o haydock_herm_resolv.o \
#	     act_herm.o act_1e.o smbubble.o smkxc.o ldriver.o adump.o invdrv.o cinvert.o coulmult.o aclip.o bfnorm.o \
#	     getband.o wvfdot.o ctrdfft.o elsdch.o jlmfft.o gradorth.o \
#	     inner.o soldrv.o subspace.o velmul.o velmuls.o lpart1.o lpart2.o lpart3.o 

SMNEWMA_OBJS = smnewma.f90 lxb_getval.o formx.o kmapr.o getwgt.o getumat.o bdallck.o uphase.o haydock.o act.o \
	     act_1e.o smbubble.o smkxc.o smlsauber.o adump.o invdrv.o cinvert.o coulmult.o aclip.o bfnorm.o \
	     getband.o wvfdot.o ctrdfft.o elsdch.o jlmfft.o gradorth.o \
	     inner.o soldrv.o subspace.o velmul.o velmuls.o  

SMNEWMA2_OBJS = lxb_getval.o formx.o kmapr.o getwgt.o getumat.o bdallck.o uphase.o haydock.o act.o \
	     act_1e.o smbubble.o smkxc.o smlsauber.o adump.o invdrv.o cinvert.o coulmult.o aclip.o bfnorm.o \
	     getband.o wvfdot.o ctrdfft.o elsdch.o jlmfft.o gradorth.o \
	     inner.o soldrv.o subspace.o velmul.o velmuls.o  

# Libraries
LIBS1 = -llapack -lblas
FFTI = -I/usr/include/
FFTWL = -L/usr/lib/x86_64-linux-gnu/ -lfftw3

DATAGEN_LIB = /home/ishiakasem/stdlib/build/gfortran_6C1C75E19A159722/stdlib/libstdlib.a
DATAGEN_INCL = -I/home/ishiakasem/stdlib/build/gfortran_020AF3609BFD7C6D

#all: $(FILES) wfconvert.x nelectron.x rhoofg.x datagen.x
all: $(PREPFILES) $(BSEFILES)

#install:
#	cp $(SCRIPTS) $(FILES) $(FILESF90) $(INSTDIR)

#instdev:
#	for F in $(SCRIPTS); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
#	for F in $(FILES); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;
#	for F in $(FILESF90); do ln -fs $(PWD)/$$F $(INSTDEVDIR); done;

.f.x:
	$(FC) $(OPTIONS) -o $*.x $*.f

# Compiling executables
#
# ----- PREP STAGE -----
tmelnew.x: tmelnew.f
	$(FC) $(OPTIONS) tmelnew.f -o tmelnew.x

whom0.x: whom0.f
	$(FC) $(OPTIONS) whom0.f -o whom0.x

frmtwf.x: frmtwf.f
	$(FC) $(OPTIONS) frmtwf.f -o frmtwf.x

rhoofg.x: rho2.f90
	$(FC) $(OPTIONS) rho2.f90 -o rhoofg.x $(FFTI) $(FFTWL)

nelectron.x: nelectron.f90
	$(FC) $(OPTIONS) nelectron.f90 -o nelectron.x

#datagen.x: $(DATAGEN_OBJS)
#	$(FC) $(OPTIONS) $(DATAGEN_LIB) $(DATAGEN_OBJS) -o datagen.x 

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

smnewma_nherm.x: smnewma_nherm.f90 $(SMNEWMA2_OBJS)
	$(FC) $(OPTIONS) smnewma_nherm.f90 $(SMNEWMA2_OBJS) -o smnewma_nherm.x $(LIBS1)

#smnewma_herm_no.x: smnewma_herm_no.f90 $(BSE_OBJS)
#	$(FC) $(OPTIONS) smnewma_herm_no.f90 $(BSE_OBJS) -o smnewma_herm_no.x $(LIBS1)

#smnewma_herm_resolv.x: smnewma_herm_resolv.f90 $(BSE_OBJS)
#	$(FC) $(OPTIONS) smnewma_herm_resolv.f90 $(BSE_OBJS) -o smnewma_herm_resolv.x $(LIBS1)

spectrum.x: spectrum.o lxb_getval.o
	$(FC) $(OPTIONS) spectrum.f lxb_getval.o -o spectrum.x

zspectrum.x: zspectrum.o lxb_getval.o
	$(FC) $(OPTIONS) zspectrum.f lxb_getval.o -o zspectrum.x

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
	$(FC) $(OPTIONS) lxb_getval.f90 -c lxb_getval.o 

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
	$(FC) $(OPTIONS) jlmfft.f -c jlmfft.o 

formx.o: formx.f90
	$(FC) $(OPTIONS) formx.f90 -c formx.o 

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
#
kmapr.o: kmapr.f90
	$(FC) $(OPTIONS) kmapr.f90 -c kmapr.o 

getwgt.o: getwgt.f90
	$(FC) $(OPTIONS) getwgt.f90 -c getwgt.o 

getumat.o: getumat.f90
	$(FC) $(OPTIONS) getumat.f90 -c getumat.o 

bdallck.o: bdallck.f90
	$(FC) $(OPTIONS) bdallck.f90 -c bdallck.o 

uphase.o: uphase.f90
	$(FC) $(OPTIONS) uphase.f90 -c uphase.o 

haydock.o: haydock.f90
	$(FC) $(OPTIONS) haydock.f90 -c haydock.o 

#haydock_herm.o: haydock_herm.f90
#	$(FC) $(OPTIONS) haydock_herm.f90 -c haydock_herm.o 

#haydock_herm_resolv.o: haydock_herm_resolv.f90
#	$(FC) $(OPTIONS) haydock_herm_resolv.f90 -c haydock_herm_resolv.o 

sizereport.o: sizereport.f90
	$(FC) $(OPTIONS) sizereport.f90 -c sizereport.o 

act.o: act.f90
	$(FC) $(OPTIONS) act.f90 -c act.o 

#act_herm.o: act_herm.f90
#	$(FC) $(OPTIONS) act_herm.f90 -c act_herm.o 

act_1e.o: act_1e.f90
	$(FC) $(OPTIONS) act_1e.f90 -c act_1e.o 

smbubble.o: smbubble.f90
	$(FC) $(OPTIONS) smbubble.f90 -c smbubble.o 

smkxc.o: smkxc.f90
	$(FC) $(OPTIONS) smkxc.f90 -c smkxc.o 

smlsauber.o: smlsauber.f90
	$(FC) $(OPTIONS) smlsauber.f90 -c smlsauber.o 

adump.o: adump.f90
	$(FC) $(OPTIONS) adump.f90 -c adump.o 

invdrv.o: invdrv.f90
	$(FC) $(OPTIONS) invdrv.f90 -c invdrv.o 

cinvert.o: cinvert.f
	$(FC) $(OPTIONS) cinvert.f -c cinvert.o 

coulmult.o: coulmult.f90
	$(FC) $(OPTIONS) coulmult.f90 -c coulmult.o 

aclip.o: aclip.f90
	$(FC) $(OPTIONS) aclip.f90 -c aclip.o 

bfnorm.o: bfnorm.f90
	$(FC) $(OPTIONS) bfnorm.f90 -c bfnorm.o 

getband.o: getband.f90
	$(FC) $(OPTIONS) getband.f90 -c getband.o 

#ldriver.o: ldriver.f90
#	$(FC) $(OPTIONS) ldriver.f90 -c ldriver.o 

lpart1.o: lpart1.f90
	$(FC) $(OPTIONS) lpart1.f90 -c lpart1.o 

lpart2.o: lpart2.f90
	$(FC) $(OPTIONS) lpart2.f90 -c lpart2.o 

lpart3.o: lpart3.f90
	$(FC) $(OPTIONS) lpart3.f90 -c lpart3.o 

tid.o: tid.f90
	$(FC) $(OPTIONS) tid.f90 -c tid.o 

wvfdot.o: wvfdot.f90
	$(FC) $(OPTIONS) wvfdot.f90 -c wvfdot.o 

ctrdfft.o: ctrdfft.f90
	$(FC) $(OPTIONS) ctrdfft.f90 -c ctrdfft.o 

elsdch.o: elsdch.f
	$(FC) $(OPTIONS) elsdch.f -c elsdch.o 

gradorth.o: gradorth.f
	$(FC) $(OPTIONS) gradorth.f -c gradorth.o 

inner.o: inner.f
	$(FC) $(OPTIONS) inner.f -c inner.o 

soldrv.o: soldrv.f
	$(FC) $(OPTIONS) soldrv.f -c soldrv.o 

subspace.o: subspace.f
	$(FC) $(OPTIONS) subspace.f -c subspace.o 

velmul.o: velmul.f90
	$(FC) $(OPTIONS) velmul.f90 -c velmul.o 

velmuls.o: velmuls.f90
	$(FC) $(OPTIONS) velmuls.f90 -c velmuls.o 

spectrum.o: spectrum.f
	$(FC) $(OPTIONS) spectrum.f -c spectrum.o 

zspectrum.o: zspectrum.f
	$(FC) $(OPTIONS) zspectrum.f -c zspectrum.o 
clean:
	rm -f *.x *.o
