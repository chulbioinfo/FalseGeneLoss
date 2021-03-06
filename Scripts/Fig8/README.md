Section for the distinguishment erroneous missing from true individual deletions.
- Zip.files of each species for data and python scripts (https://drive.google.com/drive/folders/12QbKg1SB6A_E5kz4hsv5Aw0VStAQPrNN?usp=sharing)

Accesible links for prior raw reads
- Platypus : Prior Sanger reads (https://www.ncbi.nlm.nih.gov/sra?LinkName=bioproject_sra_all&from_uid=12885)
- Climbing perch : Prior Illumina reads (https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=ERR1473920) 

Softwares
1. Read mapping: **minimap2** (v2.22-r1105-dirty) with the options: -ax map-pb for Sanger reads and -ax for Illumina reads
   (https://github.com/lh3/minimap2)
2. Read counts and visualizations: **IGV tools** 
   (https://software.broadinstitute.org/software/igv/)
3. Mppaing to assembly gaps in previous assemblies: **HAL**
   (https://github.com/ComparativeGenomicsToolkit/hal)
4. Calculating proportions of errorneous sites in missing/existing regions: **Python scripts** in each zip file
