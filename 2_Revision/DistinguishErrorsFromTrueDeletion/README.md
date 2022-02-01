Data and scripts to distinguish erroneous missing from true individual deletions.

- Species: Platypus, Climbing perch
- Prior sequening technologies: Sanger and Illumina platforms
- Zip.files of each species for data and python scripts (https://drive.google.com/drive/folders/12QbKg1SB6A_E5kz4hsv5Aw0VStAQPrNN?usp=sharing)

Softwares
1. Read mapping: **minimap2** (v2.22-r1105-dirty) with the options: -ax map-pb for Sanger reads and -ax for Illumina reads
   (https://github.com/lh3/minimap2)
2. Read counts and visualizations: **IGV tools** 
   (https://software.broadinstitute.org/software/igv/)
3. Calculating proportions of errorneous sites in missing/existing regions: **Python script**
