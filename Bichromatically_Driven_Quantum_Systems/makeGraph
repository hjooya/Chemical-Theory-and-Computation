FC = /.../intel/oneapi/compiler/2022.2.0/linux/bin/intel64/ifort
OP = -i8
LIB = -L /.../intel/oneapi/mkl/2022.2.0/lib/intel64/ \
      -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl

Bichromatically_Driven_TLSCode : Bichromatically_Driven_TLS.o
	$(FC) $(OP) -o Bichromatically_Driven_TLSCode Bichromatically_Driven_TLS.o $(LIB)
Bichromatically_Driven_TLS.o : Bichromatically_Driven_TLS.f90
	$(FC) $(OP) -c Bichromatically_Driven_TLS.f90


