cd ../scripts
inputReadDataPath = '../data/ReadData/AllDat_Chr10WT.mat';
inferredRateSavingPath = '../data/InferedRates/kineticRateChr10.mat';
[Sites, inferedRates, inferredMethyFrac] = kineticRatesMLE(inputReadDataPath,inferredRateSavingPath);