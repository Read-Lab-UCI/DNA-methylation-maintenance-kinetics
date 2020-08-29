# DNA-methylation-maintenance-kinetics

### !! IMPORTANT UPDATE
We corrected a data issue caused by site indexing problem in the MLEInference code. If you downloaded the data before **August 29, 2020**, please re-download the code and data for the corrected version.

## Abstract
This is a code repository for [Busto-Moner, et al., BioRxiv 2019](https://doi.org/10.1101/677013), which studies the DNA methylation maintenance by statistical analysis and stochastic modeling.

## Files Description

#### MLEInference.m
The matlab code for parameter inference in this paper.

#### data
ReadData: The input read data for kinetic rate inference in matlab format(.mat), which contains 2 variables: 
        
        sites: NSites x 1 integeters which represent the loci on the chromosome; 
        AllDat: this is an array of size (NSites, NTimepoints, 2). AllDat(i, j, 1) is the number of methylated reads at site i at timepoint j. AllDat(i, j, 2) is the number of unmethylated reads at site i at timepoint j.
    
InferedRates: The estimated parameters in tab separated text(BED format). Each file has the following fields:

		column number: content, meaning
		1st: chromosome, the name of the chromosome
		2nd: start, the starting position of the feature in the chromosome 
		3rd: end, the ending position of the feature in the chromosome 
		4th: MLE_rate, the estimated kinetic rates of remethylation
		5th: LB_CI_rate, the lower bound of the confident interval for the estimated rate
		6th: UB_CI_rate, the upper bound of the confident interval for the estimated rate
		7th: Edge_flag, a flag that is defined as follows: -1: sites where k is undefined, set k=0 (generally because no methylated reads-->f=0); 0: normal sites - k and full CI95 are identifiable; 1: sites on edge: only lower bound of k is identifiable. Reported k is an estimate of lower bound, corresponding to lower CI75; 2: sites where maxLL is not on boundary, but CI95 is. So k is identifiable, but lower CI95 bound is not. Reported upper CI95 is a lower-bound estimate
		8th: MLE_frac, the estimated steady state methylation fraction
		9th: LB_CI_frac, the lower bound of the confident interval for the estimated methylation fraction
		10th: UB_CI_frac, the upper bound of the confident interval for the estimated methylation fraction
		11th: NumRead_t0, the number of reads at 0h for the data in Charlton and J., Downing et al, 2018.
		12th: NumRead_later, the sum of read number for 1h, 4h and 16h.

InferedRates: The estimated parameters in matlab format(.mat) for each chromosome, and all files have 2 variables when loaded: 

		fittedSites: the loci of fitted parameter on the chromosome (N x 1, N is the number of fitted sites)
		inferedRates: the fitted kinetic rates matrix with the estimate and bounds. N * 4 matrix, 1st, 2nd, 3rd, 4th col: fitted MLE_rate, LB_CI_rate, UB_CI_rate, Edge_flag (see above)
		inferredMethyFrac: the fitted steady state methylation level and bounds. N * 3 matrix, 1st, 2nd, 3rd col: fitted MLE_frac, LB_CI_frac, UB_CI_frac.

#### methylation-kinetic-simulation
The models and input data for methylation kinetic simulation in this paper, i.e. distributive, processive, collaborative models.
