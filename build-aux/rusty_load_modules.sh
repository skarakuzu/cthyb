#%Module1.0

#set     name        devenv
#set     version     clang-py3-mkl
#set	opt	    /mnt/home/wentzell/opt

#set     description "Development environment based on clang, python3 and mkl"

#module-whatis   "$description"

#proc ModulesHelp { } {
#    global description url version git_hash
#    puts stderr "Description: $description"
#    puts stderr "URL:         $url"
#    puts stderr "Version:     $version"
#}

# Only one version of llvm can be loaded at a time
#conflict $name
#export MODULEPATH=$PWD:$MODULEPATH

module purge 
#module load modules-spack
module load tmux git vim
module load gcc/10.2.0 llvm/12.0.1
module load julia cuda
module load cmake gdb valgrind gperftools likwid swig singularity
module load intel-oneapi-mkl/2021.4.0
module load gmp gsl fftw nfft eigen
module load boost/1.77.0
module load magma
module load python/3.8.12 
module load openmpi/4.0.6
module load python-mpi/3.8.12-mpi
#module load gcc/11.2.0 llvm/11.1.0
module load hdf5/1.12.1

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

# ----
#

#setenv CC 	clang
#setenv CXX 	clang++
#setenv CXXFLAGS "-stdlib=libc++ -Wno-register -march=broadwell"
#setenv CXXFLAGS "-Wno-register -march=broadwell --gcc-toolchain=/cm/shared/sw/nix/store/xwhr9c54p7grmakfvyn5z3hiv7jwqy8x-gcc-10.2.0"
#setenv FC 	gfortran

#setenv ASAN_OPTIONS 	symbolize=1:detect_leaks=0
#setenv UBSAN_OPTIONS 	symbolize=1:print_stacktrace=1:halt_on_error=1
#setenv TSAN_OPTIONS 	symbolize=1:halt_on_error=1
#setenv MSAN_OPTIONS 	symbolize=1:halt_on_error=1
