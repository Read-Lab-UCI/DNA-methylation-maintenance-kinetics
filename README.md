### !! UPDATE
We corrected a data issue on this repo caused by site indexing problem in the MLEInference code. If you downloaded the data before **October 10, 2020**, please re-download the code and data for the corrected version. Please note that **this bug did not affect any results in the accompanying article**. The bug was introduced, and then fixed, after the article was submitted for publication. 

# DNA-methylation-maintenance-kinetics

## Abstract
This is a code repository for [Busto-Moner, et al., PLOS Comp Bio](https://doi.org/10.1371/journal.pcbi.1007195) and [Busto-Moner, et al., BioRxiv 2019](https://doi.org/10.1101/677013), which studies the DNA methylation maintenance by statistical analysis and stochastic modeling.

## Files Description

#### MLEInference.m
The matlab code for parameter inference in this paper.

#### data
ReadData: The input read data for kinetic rate inference in matlab format(.mat), which contains 2 variables: 
        
	sites: NSites x 1 integeters which represent the loci on the chromosome; 
	AllDat: this is an array of size (NSites, NTimepoints, 2). AllDat(i, j, 1) is the number of methylated reads at site i at post-replication timepoint j. AllDat(i, j, 2) is the number of unmethylated reads at site i at timepoint j. Data were obtained from Charlton, et al. 2018 (https://www.nature.com/articles/s41594-018-0046-4).
    
MLE_data: The estimated parameters in tab separated text(BED format). Each file has the following fields:

	column number: content, meaning
	1st: chromosome, the name of the chromosome
	2nd: start, the starting position of the feature in the chromosome 
	3rd: end, the ending position of the feature in the chromosome 
	4th: MLE_rate, the estimated kinetic rate of remethylation post-replication
	5th: LB_CI_rate, the lower bound of the confident interval for the estimated rate
	6th: UB_CI_rate, the upper bound of the confident interval for the estimated rate
	7th: MLE_frac, the estimated steady state methylation fraction (steady-state fraction of cells in population mehylated at that site)
	8th: LB_CI_frac, the lower bound of the confident interval for the estimated methylation fraction
	9th: UB_CI_frac, the upper bound of the confident interval for the estimated methylation fraction
	10th: NumRead_t0, the number of reads at 0h for the data in Charlton and J., Downing et al, 2018.
	11th: NumRead_later, the sum of read number for 1h, 4h and 16h.

InferedRates: The estimated parameters in matlab format(.mat) for each chromosome, and all files have 2 variables when loaded: 

	fittedSites: the loci of fitted parameter on the chromosome (N x 1, N is the number of fitted sites)
	inferedRates: the fitted kinetic rates matrix with the estimate and bounds. N * 4 matrix, 1st, 2nd, 3rd, 4th col: fitted MLE_rate, LB_CI_rate, UB_CI_rate, Edge_flag (see above)
	inferredMethyFrac: the fitted steady state methylation level and bounds. N * 3 matrix, 1st, 2nd, 3rd col: fitted MLE_frac, LB_CI_frac, UB_CI_frac.

#### methylation-kinetic-simulation
The models and input data for methylation kinetic simulation in this paper, i.e. distributive, processive, collaborative models.
