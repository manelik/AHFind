



all:	ahfind

ahfind:	AHaxis.f90
	@echo "Compiling"
	gfortran AHaxis.f90 -o ahfind

