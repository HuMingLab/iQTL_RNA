# iQTL pipline for RNA-seq
## Allelic RNA-seq mapping for CCGG.

Find the Hi-C pipeline [here](https://github.com/Xieeeee/AlleliC/).
Find the ATAC-seq pipeline [here](https://github.com/lindsayhrlee/iQTL_ATAC).

### Requirements

- STAR v2.7.9a
- Bedtools v2.29.1
- Samtools v1.9
- Python v3.6.8
    * Numpy v1.19.5
    * Pandas v1.1.5


### Running the pipeline
This pipeline can be run simply by:

```
./run_rna.sh [Mouse] [CC1] [CC2] [Name] [fastq R1] [fastq R2]
```

 The required input variables are:
1.	Mouse = The mouse ID
2.	CC1 = Maternal Genome Name
3.	CC2 = paternal Genome Name 
4.	Name = Mouse
5.	fastq R1 = FASTQ file for R1
6.	fastq R2 = FASTQ file for R2




### Contact Us

For any questions regarding this software, contact Ming Hu (hum@ccf.org) or Lindsay Lee (leeh7@ccf.org).
