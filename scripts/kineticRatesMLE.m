function [fittedSites, inferedRates, inferredMethyFrac] = kineticRatesMLE(inputReadDataPath,inferredRateSavingPath)
% kineticRatesMLE - Maximum likelihood estimation (MLE) for site-specific
% kinetic rates of DNA methylation maintenance.
% 
% Inputs:
%    inputReadDataPath - Path of input read data in .mat format
%    inferredRateSavingPath - Path for saving the inferred kinetic rates
%
% Outputs:
%    fittedSites - The location of the fitted CpG sites on genome: N integers, N
%    is the number of CpGs
% 
%    inferedRates - The inferred kinetic rates with upper and lower 
%                   confidence intervals(CIs), N by 3 matrix
%
%    inferredMethyFrac - The inferred steady state methylation fraction with upper and lower 
%                   confidence intervals(CIs), N by 3 matrix
%
% Example: 
%    kineticRatesMLE('../data/ReadData/AllDat_Chr1WT.mat', 
%                    '../data/InferedRates/kineticRateChr1.mat')

% load read data
load(inputReadDataPath, 'AllDat', 'sites'); 

% AllDat is an array of size (NSites, NTimepoints, 2).
% AllDat(i, j, 1) is the number of methylated reads at site i at timepoint j.
% AllDat(i, j, 2) is the number of unmethylated reads at site i at timepoint j.
% Totalreads(i, j) = AllDat(i, j, 1) + AllDat(i, j, 2).
% sites is an CpG location array with size (NSites, 1)

% timePoints : the experimental timepoints, and they are shifted by tShift, 
% which accounts for the fact that at, e.g., time=0 hrs, the captured reads
% started replicating within the previous hour. This is done also because 
% the model assumes zero probability of methylation occuring by time=0.

tShift = 0.5;
timePoints = [0, 1, 4, 16] + tShift;
numSites = numel(sites);

methyFracGrid = 0 : 0.01 : 1;    %grid of f-values at which LogLikelihood will be computed
rateGrid = 10.^(-2 : .01 : 1);   %grid of k-values
maxRate = max(rateGrid);
minRate = min(rateGrid);

%initialize various arrays
inferedRates = zeros(numSites, 3);
inferredMethyFrac = zeros(numSites, 3);
fittedSites = sites;

tic

%loop over the sites for rate and steady state methylation level MLE
parfor ii = 1 : numSites  %par
    %get the read data for site ii
    methylatedReadTimeCourseForSiteii = AllDat(ii, : , 1);
    unmethylatedReadTimeCourseForSiteii = AllDat(ii, : , 2);
    [inferedRates(ii, : ), inferredMethyFrac(ii, : ), LogLikelihood] = siteMLE(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii)
end
save(inferredRateSavingPath, 'fittedSites', 'inferedRates', 'inferredMethyFrac')
toc
end