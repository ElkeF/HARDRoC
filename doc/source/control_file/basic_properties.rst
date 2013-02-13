

Control File
============

In this file you specify what to calculate. The keywords are not requiered to
be present. However, you might get nonsense, if you do not specify the atom types
of your initial and final states.



Decay processes
---------------

ICD
~~~
For the calculation of ICD decay widths, the widths of all possible
pairs in the clusters needs to be evaluated::

  pairs  T

This will set the variable do_pairs to TRUE. (Default = FALSE)


ETMD3
~~~~~
For the calculation of ETMD3 decay widths, the widths of all possible
triples in the cluster need to be evaluated::

  triples T

This will activate do_triples. (Default = FALSE)



Hole specification
------------------

ICD and ETMD processes are characterized by the two letter atom types of initial and
final states. Only pairs and triples from the *coordinate file* consisting of
the atom types specified here, will be evaluated::

  atin   Ar
  atfin1 Xe
  atfin2 Xe

Hereby, **atin** is the atom type of the initially ionized atom, **atfin1** is
the atom type of the atom filling the vacancy in the inner valence and **atfin2**
is the atom type of the atom emitting the ICD or ETMD electron.

For ICD atfin1 does not need to be specified, since the hole filling from the
same atom is implied.



Example
-------

A typical *control file* looks like::

  pairs   T
  triples F
  atin    Ar
  atfin1  Xe
  atfin2  Xe

You also find an example in the folder test_inputs.

Since the programme is based on experimental values, some atom types
(at the moment almost all) have no entries in the corresponding library
files and are therefore currently not doable. Please contact me in case and
be so kind to add references for radiative lifetimes and ionzation cross section
tables. 
