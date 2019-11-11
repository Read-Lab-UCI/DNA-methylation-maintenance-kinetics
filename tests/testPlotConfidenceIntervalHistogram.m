clear;
cd ../code_archive
inferredRatePath = '../data/InferedRates/kineticRateChr1.mat';
figSavingPath = '../figures/';
maxRateCI = 2.5;
maxFracCI = 1;
plotConfidenceIntervalHistogram(inferredRatePath, figSavingPath, maxRateCI, maxFracCI)