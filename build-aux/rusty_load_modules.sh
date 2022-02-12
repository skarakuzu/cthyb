
module purge 
#module load modules-spack
module load tmux git vim
module load gcc/10.3.0
#module load julia 
#module load julia cuda
module load cmake gdb valgrind gperftools likwid swig singularity
module load intel-oneapi-mkl/2022.0.1
module load gmp gsl fftw nfft eigen
#module load boost/1.78.0
module load boost/1.78.0-libcpp
#module load magma
module load python/3.8.12 
module load openmpi/4.0.7
module load python-mpi/3.8.12-mpi
#module load gcc/11.2.0 llvm/11.1.0
module load hdf5/1.12.1
module load llvm/13.0.0

#export CC=mpicc
#export CXX=mpicxx
export CC=clang
export CXX=clang++
#export CC=/cm/shared/sw/nix/store/czwcs79bs7z5hdqkfjz50s8rxzh2h97m-llvm-11.1.0/bin/clang
#export CXX=/cm/shared/sw/nix/store/czwcs79bs7z5hdqkfjz50s8rxzh2h97m-llvm-11.1.0/bin/clang++
#export CXXFLAGS="-Wno-register -march=broadwell --gcc-toolchain=/cm/shared/sw/nix/store/7rg0h4s7yn66l9bkyzv2bgn9zbrvpjgm-gcc-10.2.0/bin/gcc"
#export CXXFLAGS="-Wno-register -march=broadwell --gcc-toolchain=/cm/shared/sw/nix/store/xwhr9c54p7grmakfvyn5z3hiv7jwqy8x-gcc-10.2.0"

export CXXFLAGS="-stdlib=libc++"
#export CXXFLAGS="-stdlib=libc++ -Wno-register -march=broadwell"

