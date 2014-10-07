FC       = gfortran
#FFLAGS   = -Og -p -fcheck=all -fbackslash -Wall -Wextra -fbounds-check -fbacktrace
FFLAGS   = -O2 -fbackslash

example.x: ../surface_physics.f90 example.f90
	$(FC) $(FFLAGS) -o $@ $^
	

run_example: example.x
	./example.x example.namelist
        
clean:
	rm -rf *.mod *.o *.x *.out
