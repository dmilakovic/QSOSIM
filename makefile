# Makefile for qsosim8
# JKW, Dec 2013
# Compiler:
FC= gfortran

# List of executables to be built within the package

# this is the default (ie 'make' will make for Linux as below)
#all: $(PROGRAMS)
all: q1 q9

obj1=	spline.f qsosim9.f readfits.f absdist.f power_laws.f generate.f
obj2=	dsepvar.o ewred.o spvoigt.o voigt.o vp_lycont.o f13_read.o \
		readfits.o qsosim9.o writefits.o spline.o absdist.o \
		power_laws.o

#  Macbook Pro version
q1: generate.f spline.f absdist.f qsosim9.f power_laws.f
	$(LINK.f) -c $(obj1)
q9: qsosim9.o $(obj)
	$(LINK.f) -o generate generate.o $(obj2) \
	-L/opt/local/lib -lpgplot -L/usr/X11/lib -lX11 \
	-L/opt/local/lib -lcfitsio	
