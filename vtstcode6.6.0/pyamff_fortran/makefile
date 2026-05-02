include makefile.include

SRC=fptype.f90 atoms.f90 neighborlist.f90 behlerParrinello.f90 nntype.f90 trainType.f90 normalization.f90 mtfort90.f90 fingerprints.f90 neuralnetwork.f90 nnmodule.f90 pyamffType.f90 lossgrad.f90 adam.f90 lbfgs_ml.f90 rprop.f90 opts.f90 training.f90 model.f90 

.PHONY: all

all: $(SRC)
	$(FC) -c $+
	cp *.mod ..
clean: 
	rm -rf *.mod *.o
