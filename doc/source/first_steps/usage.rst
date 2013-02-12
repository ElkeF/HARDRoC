

General usage
=============

HARDRoC requires four input files::

  $ ./hardroc.x control.in coordinates.xyz icd_channels etmd_channels


**Control file**
  The control file contains information about which kind of decay process
  to calculate and which atoms to choose for initial and final states.

**Coordinate file**
  The coordinate file contains an xyz file of the cluster.

**File with ICD channel specification**
  In this file all ICD channels to be investigated are listed with their
  parameters.

**File with ETMD channel specifiction**
  In this file all ETMD channels to be investigated are listed together
  with their parameters.


**Important**: All four files have to be given to the executable. Even
if no ICD calculation is performed, some file has to be given to the programme.
This file is allowed to contain anything and won't be processed if not
specified in the *control file*, but has to be present.
