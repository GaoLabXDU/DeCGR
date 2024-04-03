# DeCGR
DeCGR, an interactive tool for deciphering CGRs from chromatin contact maps. DeCGR can automatically reconstruct CGRs using the input rearrangement breakpoints and the Hi-C map. DeCGR can also validate CGRs by simulating the Hi-C map with CGRs and contrasting the abnormal interaction regions with the original contact map. Furthermore, DeCGR presents the reconstructed Hi-C contact maps, allowing visualization of the impact of rearrangements on chromatin structure after the CGRs are determined.
# Required Packages
The following software must be installed on your machine:
python (required 3.7+ and tested with version 3.9)
## Python dependencies
* iced (tested 0.5.13)
* matplotlib (tested 3.8.2)
* numpy (tested 1.26.3)
* pandas (tested 2.2.0)
* PyQt5 (tested 5.15.10) 
* scikit-learn (tested 1.4.0)
* cooler (tested 0.9.3)
* pyBigWig (tested 0.3.22)   
* sortedcontainers (tested 2.4.0)
* pyBigWig (tested 0.3.22)
# Installation
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels r
conda create -n DeCGR python=3.9 rpy2 r-mgcv
conda activate DeCGR
pip install DeCGR -U
```
# Tutorial
DeCGR consists of three modules for deciphering complex genome rearrangements from chromatin contact maps: (1) Fragment assembly, (2) Validation CGRs by simulation, and (3) Hi-C map reconstruction. Test data is provided in the TestData folder.
## Use graphical interface
$ python DeCGR.py
## Fragment assembly
### input files
* Tumor Hi-C (Hi-C file in mcool format)   
* Structural variations (The path of rearrangement breakpoints file)
For example
```
$ head K562_chr18_result.txt
chr18	450000	chr18	21850000
chr18	600000	chr18	3450000
chr18	3450000	chr18	10850000
chr18	7350000	chr18	25900000
chr18	8050000	chr18	23700000
```
* Resolution (The resolution of Hi-C matrix)
* chromosome (Chromosome ID)
### Output files
* Assembly result
```
$ head assembly_result.txt
num	chrom	start	end	node	orient
1	chr18	150000	450000	A	+
1	chr18	21850000	23700000	B	+
1	chr18	8050000	11000000	F	+
1	chr18	1650000	3450000	E	-
2	chr18	600000	850000	C	-
2	chr18	3450000	3850000	D	+
3	chr18	7350000	7950000	G	-
3	chr18	25900000	30250000	H	+
```
## Simulation
### input files
* Control Hi-C (noraml Hi-C file in mcool format)   
Tumor Hi-C (tumor Hi-C file in mcool format)
* Assembly result (All reconstructed CGRs. It can use the output of the fragment assembly as input, or the user can modify the file as input.)
* Using Fithic
When the user uses this module for the first time, the user needs to click the fithic button to run fithic. If the user runs the same data for the second time, they can use the generated dis2cnt file as input, which is stored in (../DeCGR/result/dis2cnt).
## Reconstruction
### input files
* If the user chooses cnv normalization, the user needs to provide reference and RE information. If the user chooses ice normalization, these are not needed.
* Assembly complex rearrangement (A complete complex rearrangement)
```
$ head CGR_1.txt
chrom	start	end	node	orient
chr18	150000	450000	A	+
chr18	21850000	23650000	B	+
chr18	8050000	10900000	F	+
chr18	1650000	3450000	E	-
```
