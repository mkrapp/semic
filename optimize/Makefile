.PHONY: clean run_example
FC       = gfortran
#FFLAGS   = -Og -p -fcheck=all -fbackslash -Wall -Wextra -fbounds-check -fbacktrace
FFLAGS   = -O3 -fbackslash -fcheck=all -cpp
LDFLAGS  = 
INCLUDES =
LFLAGS   =

ifeq ($(mpi),1)
	FC = mpif90
	FFLAGS += -DMPI
endif

run_particles.x: ../surface_physics.f90 ../utils.f90 run_particles.f90
	$(FC) $(FFLAGS) -o $@ $^

python-install:
	python setup.py install

clean:
	rm -rf *.mod *.o *.x *.out
