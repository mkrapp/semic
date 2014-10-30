.PHONY: clean run_example
FC       = gfortran
#FFLAGS   = -Og -p -fcheck=all -fbackslash -Wall -Wextra -fbounds-check -fbacktrace
FFLAGS   = -O2 -fbackslash -fcheck=all

example.x: ../surface_physics.f90 ../utils.f90 example.f90
	$(FC) $(FFLAGS) -o $@ $^
	
output = $(shell grep vali example.namelist | awk -F '"' '{print $$2}')
run_example: example.x
	./example.x example.namelist
	python plot_example.py example.out $(output) 0
        
clean:
	rm -rf *.mod *.o *.x *.out
