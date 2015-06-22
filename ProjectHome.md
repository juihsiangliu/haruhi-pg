# Description #

It is a parallel power grid simulator.

# Release #

The first version is released at 2011-04-06.
Download here - http://code.google.com/p/haruhi-pg/downloads/list

## Install ##

### With library files ###
The pre-compiled third party library files are included (in Ubuntu 8.04 environment).
```
tar zxvf haruhi-pg.tar.gz
```

### Without library files ###
```
tar zxvf haruhi-pg-wo.tar.gz
```

After un-tar the file.
You may see four directories - bin, src, case, and lib.

I did not add the libraries to here, since the size is too large and they are OS dependent.
You can download them and install them in the _lib_ directory.
For more details, you could reference the _src/Makefile_ to see the directory name setting.
  * AMD - http://www.cise.ufl.edu/research/sparse/amd/
  * metis - http://glaros.dtc.umn.edu/gkhome/views/metis
  * gsl - http://www.gnu.org/software/gsl/
  * gdsl - http://home.gna.org/gdsl/
  * ufconfig - http://www.cise.ufl.edu/research/sparse/UFconfig/

## Compile ##

Enter _src/_ and just type _make_.
If all libraries are set up in the correct path, the binary file **haruhi-pg** will be generated in _bin_ directory.

## Bin ##

A pre-compiled static binary **haruhi-pg-static** is put here.
It may executed directly in Linux environment.

## Case ##

In this directory, I put _ibmpg1_ case here.
All other cases can be found in - http://dropzone.tamu.edu/~pli/PGBench/

## Execute ##

You may first enter the case directory _case/ibmpg1_, and type
```
../../bin/haruhi-pg-static ibmpg1.spice > sol.txt
```
to run the program in default mode.

Or you may use different number of thread to solve.
Just type
```
../../bin/haruhi-pg-static ibmpg1.spice -thread 1 -order_method metis > sol.txt
../../bin/haruhi-pg-static ibmpg1.spice -thread 2 -order_method metis > sol.txt
../../bin/haruhi-pg-static ibmpg1.spice -thread 4 -order_method metis > sol.txt
```
to run with different number of thread.