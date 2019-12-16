function [fittedSites, inferedRates, inferredMethyFrac, AllDat] = kineticRatesMLE(inputReadDataPath,inferredRateSavingPath)
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
%                   confidence intervals(CIs), N by 4 matrix
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
rateGrid = 10.^(-2 : .01 : 1);   %grid of rate-values (denoted as k hereafter)

%initialize various arrays
inferedRates = zeros(numSites, 4); %The inferred kinetic rates with upper, lower confidence intervals(CIs) and Edge Flag, N by 4 matrix

%The Edge Flag(4th field) is defined as follows:
%-1: sites where k is undefined, set k=0 (generally because no methylated reads-->f=0)
%0: normal sites - k and full CI95 are identifiable
%1: sites on edge: only lower bound of k is identifiable.
%Reported k is an estimate of lower bound, corresponding to lower CI75
%2: sites where maxLL is not on boundary, but CI95 is. So k is
%identifiable, but lower CI95 bound is not. Reported upper CI95 is a lower-bound estimate 

inferredMethyFrac = zeros(numSites, 3); %The inferred steady state methylation fraction with upper and lower confidence intervals(CIs), N by 3 matrix
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