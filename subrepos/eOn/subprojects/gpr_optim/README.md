# gpr_dimer

Port of the GPR-dimer code from Matlab to C++

# Setup

## Conda / Micromamba

``` bash
micromamba create -p $(pwd)/.tmp -c conda-forge compilers meson cmake gfortran eigen==3.4.0 pkg-config gtest
micromamba activate $(pwd)/.tmp
```

Or perhaps if this is used very often, install the environment globally.

``` bash
micromamba env create -f environment.yml
micromamba activate gprd
```

Alternatively, `pixi shell` can be used.

# Meson Usage

If for some reason no envionment manager is used, then one can manually setup
`gtests` as below.

1. Setup `gtests` (once)

``` bash
mkdir subprojects
meson wrap install gtest
```

Note that optionally, if not using `nix` or something else, `meson wrap install eigen` will also help

2. Compile (`input` is copied over at configuration)

``` bash
meson setup builddir
meson compile -C builddir
```

Recall that a release build with `meson` is via:

``` bash
meson setup builddir --buildtype=release --optimization=3
```

3. Profit

``` bash
cd builddir
./gprd
# Sample output
Final convergence obtained after 2 relaxation phases (total number of image evaluations: 7).
Elapsed time: 18.092s

8.98514 9.948 7.88447 7.64819 9.94644 7.88399 


-0.000319603 


0.00123983 0.0022487 0.00325383 -1.93536e-07 -0.00136041 0.000147053 
```

# HOWTO [Older]

1. Download `googletests` framework:

```bash
$ git submodule init
$ git submodule update
```

2. Patch `CMakeLists.txt`:

```bash
$ cd googletest/googletest
$ patch < ../../cmake_gtest.patch
```

3. Build `googletests` framework (replace `g++` with a compiler of your choice):

```bash
$ mkdir -p build
$ cd build
$ cmake -DCMAKE_CXX_COMPILER="g++" -DCMAKE_CXX_FLAGS="-std=c++11" ../
$ make all
```

4. Go to the main directory of this repository and call `make`:

```bash
$ make OBJDIR=obj EXE_NAME=gpr.out LIBS_PATHS=-L./googletest/googletest/build/lib
```

5. Compile the internal Fortran code (separate compilation for now) and copy `mEAMCUH2` executable to the main directory:

```bash
$ cd gpr/eam_potential_fortran
$ make
$ cp mEAMCUH2 ../../
```

6. Run the code:

```bash
$ ./gpr.out
```

7. Requires Eigen3 library, download the latest version.

# Meson Alternative
This will download the dependencies without manual intervention.

``` bash
mkdir subprojects
meson wrap install eigen 
meson wrap install gtest
meson setup builddir
meson compile -C builddir
```


# Documentation

To build and view the documentation, we have to obtain the theme and tag files, this is done on the CI, but locally needs to be carried out by the user:

```bash
mkdir -p apidocs/tags
cd apidocs/tags
curl https://upload.cppreference.com/mwiki/images/f/f8/cppreference-doxygen-web.tag.xml -o cppreference-doxygen-web.tag.xml
cd ../ # now in doxygen
wget  "https://github.com/HaoZeke/doxyYoda/releases/download/0.0.2/doxyYoda_0.0.2.tar.gz"
tar xf doxyYoda_0.0.2.tar.gz
cd ../ # now at project root
doxygen apidocs/Doxygen-prj.cfg
```

Note that:
- We expect `doxygen` to be at **1.9.1**
  - Otherwise, leave the builds to the CI
- Use Javadoc

Once built, they can be perused with any HTTP/S server:

``` bash
# Either
darkhttpd html
# OR
python -m http.server 8888 --directory html
```

