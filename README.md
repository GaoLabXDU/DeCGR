# DeCGR
An interactive tool for deciphering CGRs from chromatin contact maps.
DeCGR, an interactive tool for deciphering CGRs from chromatin contact maps. DeCGR can automatically reconstruct CGRs using the input rearrangement breakpoints and the Hi-C map. DeCGR can also validate CGRs by simulating the Hi-C map with CGRs and contrasting the abnormal interaction regions with the original contact map. Furthermore, DeCGR presents the reconstructed Hi-C contact maps, allowing visualization of the impact of rearrangements on chromatin structure after the CGRs are determined.
# Required Packages
The following software must be installed on your machine:
python: tested with version 3.7
## Python dependencies
* iced (tested 0.5.10)
* matplotlib (tested 3.3.2)
* numpy (tested 1.21.5)
* pandas (tested 1.3.5)
* Pillow (tested 9.5.0)
* PyQt5 (tested 5.15.4) 
* scikit-learn (tested 1.0.2)
* scipy (tested 1.7.3)
* sortedcontainers (tested 2.4.0)
# Installation
```
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda create -n DeCGR python=3.7
conda activate DeCGR
pip install iced matplotlib numpy pandas Pillow PyQt5 scikit-learn scipy sortedcontainers
```
# Tutorial
DeCGR consists of three modules for deciphering complex genome rearrangements from chromatin contact maps: (1) Fragment assembly, (2) Validation CGRs by simulation, and (3) Hi-C map reconstruction. Test data is provided in the Data folder.
## Use graphical interface
$ python mainPlatform.py
## Fragment assembly
### input files
* Tumor Hi-C (The path of tumor Hi-C matrix folder, Hi-C files in the folder need to be named in the format _chrA_chrB_)  
For example:  
tumor Hi-C matrix folder: /mnt/d/DeCGR/Data/BM178"  
Hi-C files are named: BM178_100kb_chr12_chr12_matrix.txt  
* Structural variations (The path of rearrangement breakpoints file)
For example  
```
$ head BM178_chr12_result.txt
chr12	60600000	chr12	83800000
chr12	46500000	chr12	82300000
chr12	46200000	chr12	82100000
chr12	45800000	chr12	91200000
chr12	57750000	chr12	82800000
chr12	58550000	chr12	78000000
```
* Resolution (The resolution of Hi-C matrix)
* chromosome (Chromosome ID, it need to be consistent with the Hi-C matrix storage in the Hi-C folder)
### Output files
* Assembly result
```
$ head assembly_result.txt
num	chrom	start	end	node	orient
1	chr12	45800000	46100000	G	-
1	chr12	83800000	91200000	B	-
1	chr12	58300000	60800000	A	-
1	chr12	78000000	82300000	D	+
1	chr12	46500000	48600000	C	+
2	chr12	46200000	46500000	E	-
2	chr12	82100000	82800000	F	+
2	chr12	55900000	57700000	H	-

```
## Simulation
### input files
* Control Hi-C (The path of tumor Hi-C matrix folder, Hi-C files in the folder need to be named in the format _chrA_chrB_)  
For example:   
tumor Hi-C matrix folder: /mnt/d/DeCGR/Data/GM12878"  
Hi-C files are named: GM12878_100kb_chr12_chr12_matrix.txt  
* Assembly result (All reconstructed CGRs. It can use the output of the fragment assembly as input, or the user can modify the file as input.)
* Using Fithic
When the user uses this module for the first time, the user needs to click the fithic button to run fithic. If the user runs the same data for the second time, they can use the generated dis2cnt file as input, which is stored in (../DeCGR/result/dis2cnt).
## Reconstruction
### input files
* Assembly complex rearrangement (A complete complex rearrangement)
```
$ head CGR_1.txt
chrom	start	end	node	orient
chr12	45800000	46100000	G	-
chr12	83800000	91200000	B	-
chr12	58300000	60800000	A	-
chr12	78000000	82300000	D	+
chr12	46500000	48600000	C	+
```
