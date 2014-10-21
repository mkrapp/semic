.PHONY: clean run_example
FC       = gfortran
#FFLAGS   = -Og -p -fcheck=all -fbackslash -Wall -Wextra -fbounds-check -fbacktrace
FFLAGS   = -O2 -fbackslash -fcheck=all
LDFLAGS  = 
INCLUDES =
LFLAGS   =

ifeq ($(mpi),1)
	FC = mpif90
	FFLAGS += -DMPI
endif

run_particles.x: ../surface_physics.f90 ../utils.f90 run_particles.F90
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -rf *.mod *.o *.x *.out
