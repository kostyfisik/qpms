Notes on installing and running qpms on Triton
==============================================

Running anything in a cluster environment can cause some unexpected problems.
With these notes, one should be able to run qpms successfully on Aalto's Triton.

The main problem is that it is often not completely trivial to ensure 
mutual compatibility of various libraries and tools on which qpms depends.

## Triton+anaconda
### Installing qpms with anaconda

This section shows how to use qpms in an anaconda virtual environment.
First, we purge the modules to ensure that no other other conflicting python instances are loaded.
Then we load an anaconda3 module.
```
module purge
module load anaconda3
```
Next, conda needs to infect our shell with some code in order to work properly.
(If you use a different shell than bash, modify the following command accordingly.)
```
conda init bash
```
This adds several lines to `~/.bashrc` (or elsewhere, depending on the shell).
Then we might be required to restart the shell.

Next, we create and activate a virtual environment to be used with qpms.
It will be called `trqpms` throughout this document.
We will also need to have `cython` installed before even trying to install qpms.
```
conda create trqpms
conda activate trqpms
conda install cython
```

Moreover, we need to load the GCC and gsl modules.  GCC has to be loaded *after*
loading anaconda3 and activating the conda virtual environment, otherwise
anaconda will try to use its own (likely broken) compiler when installing qpms.
And we also need cmake for building amos.
Lastly, we save the list of loaded modules for later use as `qpms_conda`.
```
module load gsl
module load GCC
module load cmake
module save qpms_conda
```

From now on, we need to be in the qpms source root directory.
```
cd /path/to/qpms
```

Now it's time to build `amos`. First, we clean all the files previously
created by `cmake` (if needed). This is to ensure that `cmake` uses the same
`gfortran` that is in the current `PATH`; otherwise, `cmake` might use some
Fortran compiler detected earlier (`cmake` apparently does not update
that information) and that might later cause mismatch between
libgfortran versions.

```
rm -r CMakeFiles CMakeCache.txt
cmake .
make clean
make amos
```
Cmake builds in other directories are not supported right now, as qpms setup scripts now expects
certain amos files at fixed relative paths. Therefore the `.` in `cmake .`.

And finaly, build qpms.
```
python3 setup.py install
```

At this point, qpms should be installed successfully. If not, this document
needs updates. 

### Running qpms with anaconda

After installing the qpms python module, one should be able to import it and to use
the scripts in the misc directory in the same shell session used for installing 
qpms.

To use qpms also later, one has to restore the environment, i.e. load the same modules and
activate the anaconda virtual environment.
```
module restore qpms_conda
conda activate trqpms
```


