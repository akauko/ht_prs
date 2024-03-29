# ht_prs
Code for the publication: Vaura, F., et al. "Polygenic risk scores predict hypertension onset and cardiovascular risk." Hypertension 77.4 (2021): 1119-1127.
https://doi.org/10.1161/hypertensionaha.120.16471

* Data: FinnGen (https://www.finngen.fi/en)
* PRS values were calculuted for FinnGen individuals using PRS-CS pipeline with default settings: https://github.com/getian107/PRScs
* We used GWAS summaries from UKBB GWAS v3: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=227859291

```
ht_prs
├── README.md                 		# Overview
├── ht_prs_final.Rmd          		# R markdown for the analysis with FinnGen data
├── ht_prs_final_1408.html     		# Documentation html generated from R markdown
├── articles-functions.R      		# Minor R functions for the main analysis
├── select_columns.pl         		# Perl script to select columns from tsv files by column name
	
```
