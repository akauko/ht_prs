# ht_prs
R code for the publication: Polygenic Risk Scores Predict Hypertension Onset and Cardiovascular Risk
https://doi.org/10.1161/hypertensionaha.120.16471

* PRS values were calculuted for FinnGen individuals using PRS-CS pipeline with default settings: https://github.com/getian107/PRScs
* We used GWAS summaries from UKBB GWAS v3: https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=227859291

```
ht_prs
├── README.md                 		# Overview
├── ht_prs_final.rmd          		# R markdown for the analysis with FinnGen data
├── articles-functions.R      		# Minor R functions for the main analysis
├── select_columns.pl         		# Perl script to select columns from tsv files by column name
	
```
