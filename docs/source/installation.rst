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
* wget (tested 3.2)

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

Notes for Installation on Windows Systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1. **Avoid PowerShell**: Use Command Prompt instead of PowerShell for the installation, as certain commands may not work correctly in PowerShell.

2. **Installation Error**: Encountering the error "error: Microsoft Visual C++ 14.0 or greater is required."

   **Solution**: Refer to this `guide <https://stackoverflow.com/questions/64261546/how-to-solve-error-microsoft-visual-c-14-0-or-greater-is-required-when-inst>`_.

Running
-------

Use the command below to access the DeCGR interface.

Mac OS or Linux (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: bash

   $ DeCGR

Windows
~~~~~~~
.. code-block:: bash

   $ DeCGR_windows

Notes for Running on Windows Systems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On Windows systems, packages like **rpy2** may not be fully compatible, which can lead to various issues during runtime. Below are common problems and their solutions for user reference:

1. **JIT Initialization Error**: Encountering "fatal error unable to initialize the jit" when executing **DeCGR_windows**.

   **Solution**:

   .. code-block:: bash

      $ set R_HOME=%CONDA_PREFIX%\Lib\R

2. **Missing 'stats' Package**: Error "package 'stats' in options('defaultPackages') was not found" when running the Fragment Assembly or visualization modules.

   **Solution**:

   .. code-block:: bash

      $ set ENV_PATH=%CONDA_PREFIX%

      $ xcopy "%ENV_PATH%\Lib\R\bin\x64\*.*" "%ENV_PATH%\Lib\R\library\stats\libs\x64" /S /I


Test Data
---------
**Tumor Hi-C data**: K562 cell line (4DNFI18UHVRO.mcool) can be downloaded from `4dnucleome <https://data.4dnucleome.org/files-processed/4DNFI18UHVRO/>`_.

**Control Hi-C data**: GM12878 cell line (4DNFITRVKRPA.mcool) can be downloaded from `4dnucleome <https://data.4dnucleome.org/files-processed/4DNFITRVKRPA/>`_.

**Breakpoint file (for Module2)**: This file is named "K562_chr18_breakpoint.txt" and provided in the `DeCGR <https://github.com/GaoLabXDU/DeCGR/tree/main/TestData>`_.

**Assembly result file (for Module3)**: This file is named "K562_chr18_assembly_result.txt" and provided in the `DeCGR <https://github.com/GaoLabXDU/DeCGR/tree/main/TestData>`_.

**Complex rearrangement file (for Module4)**: This file is named "K562_chr18_CGR.txt" and provided in the `DeCGR <https://github.com/GaoLabXDU/DeCGR/tree/main/TestData>`_.
