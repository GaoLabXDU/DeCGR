Installation
============

Required Packages
-----------------
The following software must be installed on your machine:

* Python (required 3.7+, tested with version 3.9)

Python dependencies
~~~~~~~~~~~~~~~~~~~
* matplotlib (tested 3.9.2)
* numpy (tested 1.26.4)
* pandas (tested 2.2.3)
* PyQt5 (tested 5.15.11)
* scikit-learn (tested 1.5.2)
* cooler (tested 0.10.2)
* sortedcontainers (tested 2.4.0)
* pomegranate (tested 0.15.0)
* scikit-image (tested 0.24.0)

Installation
------------
DeCGR and all dependencies can be installed using pip::

   $ conda config --add channels r
   $ conda config --add channels defaults
   $ conda config --add channels bioconda
   $ conda config --add channels conda-forge
   $ conda create -n DeCGR python=3.9 rpy2 r-mgcv
   $ conda activate DeCGR
   $ pip install -U DeCGR

Test Data
---------
**Tumor Hi-C data**: K562 cell line (4DNFI18UHVRO.mcool) can be downloaded from `4dnucleome <https://data.4dnucleome.org/files-processed/4DNFI18UHVRO/>`_.

**Control Hi-C data**: GM12878 cell line (4DNFITRVKRPA.mcool) can be downloaded from `4dnucleome <https://data.4dnucleome.org/files-processed/4DNFITRVKRPA/>`_.

**Breakpoint file (for Module2)**: This file is named "K562_chr18_breakpoint.txt" and provided in the `DeCGR <https://github.com/GaoLabXDU/DeCGR/tree/main/TestData>`_.

**Assembly result file (for Module3)**: This file is named "K562_chr18_assembly_result.txt" and provided in the `DeCGR <https://github.com/GaoLabXDU/DeCGR/tree/main/TestData>`_.

**Complex rearrangement file (for Module4)**: This file is named "K562_chr18_CGR.txt" and provided in the `DeCGR <https://github.com/GaoLabXDU/DeCGR/tree/main/TestData>`_.