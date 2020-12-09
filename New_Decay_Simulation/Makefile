FC=gfortran
FFLAGS=-I.

#ODIR=/Users/abhijitbaishya/Documents/Work/Simulation/simulation_hoyle/New_Decay_Simulation
SDIR= .
#SDIR=/home/abhijit/Documents/New_Decay_Simulation

#OBJ = mod_interfaces.o mod_functional.o func_subs.o
OBJ = func_subs.o
#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))
 
%.o: $(SDIR)/%.f90
	$(FC) -o $@ -c $< $(FFLAGS)

decay: $(OBJ)
	$(FC) -o $@ test.f90 $^ $(FFLAGS)

.PHONY : clean

clean:
	@rm -f *.o *.mod decay

