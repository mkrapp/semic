.PHONY: clean run_example
FC       = gfortran
#FFLAGS   = -Og -p -fcheck=all -fbackslash -Wall -Wextra -fbounds-check -fbacktrace
FFLAGS   = -O2 -fbackslash -fcheck=all

run_particles.x: ../surface_physics.f90 ../utils.f90 run_particles.f90
	$(FC) $(FFLAGS) -o $@ $^

clean:
	rm -rf *.mod *.o *.x *.out
