#OURS = main.f90 struc.f90 band.f90 output.f90
#BJS = main.o struc.o band.o output.o
#IN = main
#90 = ifort

#{BIN}:${OBJS}
#{F90} -o ${BIN} ${OBJS}
#{OBJS}: %.o :%.F90
#{F90} -c $<

FC = ifort

SUFFIX = .f90

MODULE = const.o struc.o edge.o utility.o fermisurf.o test2.o
EXE = test2

default : $(EXE)

$(EXE):%: $(MODULE) %.o
	$(FC) -o $@ $^ -llapack

%.o:%.f90
	$(FC) -c $*.f90

clean:
	rm -f *.o *.mod
