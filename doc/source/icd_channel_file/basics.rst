

ICD Channel File
================

The *ICD channel file* is a collection of numbers characterising the
ICD channel to be investigated. A typical file would look like::

  # J_A   M_A  SIP(in)  shift(in)   J_Ap   M_Ap  SIP(fin1)  shift(fin1)  tauabs  taurel  j_Bp  SIP(fin2)  shift(fin2)  sigma_rel
   1     1   29.2      -0.636      3      88     12.5        -1.3       4.2E-9    2.0    3    12.5        -1.3           0.625
   1     1   29.2      -0.636      1      1      12.5        -1.3       4.2E-9    2.0    3    12.5        -1.3           1.6  

It is crucial to have the correct number of entries in every row.
If you do not need them, set the values to 0 or to 1.



Initial state characteristics
-----------------------------

J_A
~~~
This is twice the total angular momentum of the initially ionized state.
In the example given above J_A = 1/2.

M_A
~~~
This is twice the projection of J_A.

SIP(in)
~~~~~~~
is the single ionization potential of the initial ionized atom in eV.

shift(in)
~~~~~~~~~
If the ionization energy is shifted due to the surroundings, the shift in eV is
to be given here. If no shift is needed, type a 0.


Final state characteristics of the initially ionized atom
---------------------------------------------------------

J_Ap
~~~~
is twice the total angular momentum of the final state on the initially
ionized atom.

M_Ap
~~~~
is twice the projection of J_Ap. Here there is one exception: If M_Ap = 88,
then for the given value of J_Ap the decay widths of all possible values of
M_Ap are calculated and added up.

SIP(fin1)
~~~~~~~~~
Single ionization potential of the initially ionized atom.

shift(fin1)
~~~~~~~~~~~
Shift of the single ionization potential in eV.

tauabs
~~~~~~
The radiative lifetime of the state. Very often this is given for more than one
J value. Therefore taurel is needed as well.

taurel
~~~~~~
In case of tauabs including the radiative lifetime of more than one J value this
parameter is the relation between the lifetime of the here specified state and
the other state.


Final state characteristics of the electron emitting atom
---------------------------------------------------------

j_Bp
~~~~
is twice the total angular momentum of the ICD electron emitting atom.

SIP(fin2)
~~~~~~~~~
Single ionization potential of the electron emitting atom in eV.

shift(fin2)
~~~~~~~~~~~
Shift of SIP(fin2) in eV.

sigma_rel
~~~~~~~~~
Like the radiative lifetime, the ionization cross section is experimentally
achieved for more than one state. The relation between the two has to be
specified in sigma_rel. However, the ratio needs to be given as others/wanted.
This unintuitive input allows calculations with only one specific ionization
cross section as in the relativistic case. Otherwise dividing by 0 would give
unintended errors.

The total value of sigma is determined from experimental
data in a library module.
