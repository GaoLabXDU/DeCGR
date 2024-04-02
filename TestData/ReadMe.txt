Description of test data

# Tumor Hi-C data 
K562 cell line (4DNFI18UHVRO.mcool) can be downloaded from https://data.4dnucleome.org/files-processed/4DNFI18UHVRO/.

# Control Hi-C data
GM12878 cell line (4DNFITRVKRPA.mcool) can be downloaded from https://data.4dnucleome.org/files-processed/4DNFITRVKRPA/.

# Breakpoint data
This file is named "K562_chr18_result.txt" and provided in the TestData.

# Fithic file
The distance-interaction file obtained after running Fithic on GM12878 cell line. This file is named "dis2cnt" and provided in the TestData.

# CNV normalization files
* hg38.MboI.npz (This file needs to be decompressed);
* hg38_1kb_GC.bw;
* hg38_mappability_100mer.1kb.bw;
These files are required for cnv normalization of the test dataset. Storing these files in the "./result" path saves time downloading these files when running your code.

# assembly result
This file is named "assembly_result.txt" and provided in the TestData.

# complex rearrangement file
This file is named "CGR_1.txt" and provided in the TestData.






