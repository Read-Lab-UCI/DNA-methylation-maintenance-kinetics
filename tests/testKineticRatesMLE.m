cd ../scripts
inputReadDataPath = '../data/ReadData/AllDat_Chr1WT.mat';
inferredRateSavingPath = '../data/InferedRates/kineticRateChr1.mat';
[Sites, inferedRates, inferredMethyFrac] = kineticRatesMLE(inputReadDataPath,inferredRateSavingPath)