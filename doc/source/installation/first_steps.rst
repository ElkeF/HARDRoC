

Getting the Code
================

The code is under version control using git. You can get the code
by cloning it into your favorite place::

  $ git clone /home/elke/pub/hardroc your_favorite_directory



Compiling
=========

HARDRoC is written in Fortran 2003. Please make sure you have a compiler,
which is able to handle it. In case you are using gfortran, you can simply
use::

  $ ./compile

Otherwise exchange the compiler in the file "compile" by your favorite
Fortran compiler.

The executable and all the compiled modules used are created in the "src"
directory.
