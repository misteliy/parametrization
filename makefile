shell   := /bin/sh
#ftn := gfortran
ftn := ftn -g -O3 -fbounds-check 
#ftn := mpif90 -lblas -llapack -lfftw3 -lblacs -lblacsf77
#ftn := mpif90 -g -O3 -fbounds-check
#ftn := mpif90 -O3 -fbounds-check
#ftn := ifort 
$(shell rm *.o)
$(shell rm *.mod)
MYFLAG := /project/s238/misteliy/cp2k/makefiles/../lib/CRAY-XE6-gfortran/popt/libcp2k_lib.a /project/s238/misteliy/cp2k/makefiles/../lib/CRAY-XE6-gfortran/popt/libcp2k_dbcsr_lib.a /project/s238/misteliy/cp2k/makefiles/../lib/CRAY-XE6-gfortran/popt/libcp2k_fft_lib.a /project/s238/misteliy/cp2k/makefiles/../lib/CRAY-XE6-gfortran/popt/libcp2k_base_lib.a /project/s238/misteliy/cp2k/makefiles/../lib/CRAY-XE6-gfortran/popt/libcp2k_ma_lib.a /project/s238/misteliy/cp2k/makefiles/../lib/CRAY-XE6-gfortran/popt/libcp2k_elpa_lib.a    -L/users/vondele/libsmm/ -lsmm_dnn_cce_interlagos  -L/users/misteliy/libint-1.1.4/lib -lderiv -lint -lstdc++
all: hint main
main: main.o random.o obj_fun.o write_nc.o ls_rmsd.o optim.o defs.o
	${ftn} -o main main.o random.o obj_fun.o write_nc.o ls_rmsd.o optim.o defs.o ${MYFLAG} 
defs.o: defs.f90
	${ftn} -c defs.f90
main.o: main.f90 obj_fun.o random.o write_nc.o optim.o defs.o
	${ftn} -c main.f90 
ls_rmsd.o: ls_rmsd.f90 defs.o
	${ftn} -c ls_rmsd.f90
obj_fun.o:obj_fun.f90 random.o ls_rmsd.o defs.o
	${ftn} -c obj_fun.f90 
random.o:random.f90 defs.o
	${ftn} -c random.f90
optim.o:optim.f90 defs.o
	${ftn} -c optim.f90
write_nc.o:write_nc.f90 optim.o defs.o
	${ftn} -c write_nc.f90 
hint:
#	echo "####on opt3 do first!!   module load gcc/4.5.2 openmpi/1.4.3-gcc_4.5.2  intel  ####"
#	echo "==================================================================================="
	
