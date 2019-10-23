clear;
cd ../code_archive
inferredRatePath = '../data/InferedRates/kineticRateChr1.mat';
figSavingPath = '../figures/';
maxRateCI = 5;
maxFracCI = 0.5;
plotConfidenceIntervalHistogram(inferredRatePath, figSavingPath, maxRateCI, maxFracCI)