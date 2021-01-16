# comopiler
FC = ifort

# compile flags
FCFLAGS = -c

# source files and objects
SRCS = $(patsubst %.f90, %.o, %.mod)

# program name
PORGRAM = aa

all: $(PORGRAM)

$(PORGRAM): $(SRCS)
	$(FC) $(FCFLAGS) $@ $<

%.o: %.f90
    $(FC) $(FLFLAGS) -o $@ $<

 %.mod
    $(FC) $(FLFLAGS) -o $@ $<

 clean:
    rm -f *.o *.mod
