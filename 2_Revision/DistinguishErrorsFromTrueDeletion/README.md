Data and scripts to distinguish erroneous missing from true individual deletions.

- Species: Platypus, Climbing perch
- Prior sequening technologies: Sanger and Illumina platforms

Softwares
1. Read mapping: **minimap2** (v2.22-r1105-dirty) with the options: -ax map-pb for Sanger reads and -ax for Illumina reads
   ()
2. Read counts and visualizations: **IGV tools** 
   (https://software.broadinstitute.org/software/igv/)
3. Calculating proportions of errorneous sites in missing/existing regions: **Python script**
