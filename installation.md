## Create a conda environment with all dependencies:
```
mamba create -n rohan -c conda-forge \
  cmake autoconf automake libtool libpng pkg-config \
  gsl bzip2 xz zlib curl ncurses \
  gcc_linux-64 gxx_linux-64 binutils_linux-64
```
## Activate the environment:
```
conda activate rohan
```
## Create symbolic links between ncurses and curses libraries
```
ln -s $CONDA_PREFIX/lib/libncurses.so $CONDA_PREFIX/lib/libcurses.so
ln -s $CONDA_PREFIX/lib/libncurses.a $CONDA_PREFIX/lib/libcurses.a
ln -s $CONDA_PREFIX/lib/libncurses.so.6 $CONDA_PREFIX/lib/libcurses.so.6
```
## Then compile:
```
cd ROHan
make
```
