



all:	ahaxis

ahaxis:	AHaxis.f90
	@echo "Compiling"
	gfortran AHaxis.f90 -o ahaxis
