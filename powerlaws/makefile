#makefile for plot.f
FC= gfortran

pgplot=-L/opt/local/lib -lpgplot
X11=-L/usr/X11/lib -lX11

obj1=	power_laws.f spline.f plot.f
obj2=	power_laws.o spline.o plot.o

all: q1 q2
q1:  power_laws.f spline.f plot.f
	$(LINK.f) -c $(obj1)
q2:  power_laws.o spline.o plot.o
	$(LINK.f) -o plot $(obj2) $(pgplot) $(X11)

