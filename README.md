# Scripts used to look for evidence of parallel selection in stickleback


## Resources

[Popoolation](https://sourceforge.net/p/popoolation/wiki/Main/)

[Popoolation2](https://sourceforge.net/p/popoolation2/wiki/Main/)

## Process and relevant files

1. **Outlier scan** : 

file = fst_outliers_detection.R (output logical data frame in wide format).

2. **Find outliers in multiple comparisons**: 

files = FST_outliers_output_example_for_overlap_scan.tsv and outlier_overlap_script

3. **Get positions of 'parallel outliers' found in step 2**: 

file = fst_outliers_detection.R (bottom)

4. **Find genic SNPs or SNPs near genes**: 

files = gene_hits.awk, genic_snps_script, stick_prot_cod_names_tab.tsv, fst_outlier_posistion_2up.tsv (output from step 3)
