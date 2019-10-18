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

logMinRate = -2; %log10 minimum allowed fit rate
logMaxRate = 1;  %log10 maximume allowed fit rate--calculated according to ReadDepth
methyFracGrid = 0 : 0.01 : 1;    %grid of f-values at which LogLikelihood will be computed
rateGrid = 10.^(logMinRate : .01 : logMaxRate);   %grid of k-values
maxRate = max(rateGrid);
minRate = min(rateGrid);

%initialize various arrays
inferedRates = zeros(numSites, 3);
inferredMethyFrac = zeros(numSites, 3);
fittedSites = sites;

tic
savingForEveryNSites = 50000;
savingCounter = 0;

%loop over the sites for rate and steady state methylation level MLE
parfor ii = 1 : numSites 
    
    %save the fitted data periodically
%     if (ii - savingCounter * savingForEveryNSites) / savingForEveryNSites >=  1
%         save ii ii
%         savingCounter = savingCounter+1;
%         save(inferredRateSavingPath, 'fittedSites', 'inferedRates', 'inferredMethyFrac')
%         toc
%     end
    
    %get the read data for site ii
    methylatedReadTimeCourseForSiteii = AllDat(ii, : , 1);
    unmethylatedReadTimeCourseForSiteii = AllDat(ii, : , 2);
    
    % probability array of the observed methylated (or unmethylated) reads at CpG ii in data given the  grid of kinetic rates and steady-state methylation levels
    pMethyArray = zeros(numel(methyFracGrid), numel(rateGrid), numel(timePoints));
    pUnmethyArray = zeros(numel(methyFracGrid), numel(rateGrid), numel(timePoints));
    
    %loop over the timepoints to compute the LogLikelihood surface
    for tind = 1 : numel(timePoints)
        timeI = timePoints(tind);
        expoentialTerms = 1 - exp(-rateGrid .* timeI);
        expoentialTerms(expoentialTerms >= 1) = 1 - eps;
        repExpoentialTerms = repmat(expoentialTerms, numel(methyFracGrid), 1);
        repMethyFrac = repmat(methyFracGrid', 1, numel(rateGrid));
        pMethyArray( : , : , tind) = methylatedReadTimeCourseForSiteii(tind) .* log(repMethyFrac .* repExpoentialTerms);
        pUnmethyArray( : , : , tind) = unmethylatedReadTimeCourseForSiteii(tind) .* log(1 - repMethyFrac .* repExpoentialTerms);
    end
    
    %this is the LogLikelihood surface as a function of the parameters
    LogLikelihood = sum(pMethyArray, 3, 'omitnan') + sum(pUnmethyArray, 3, 'omitnan');
    
    %find the parameter values that maximize the LogLikelihood
    szL = size(LogLikelihood); % size of the likelihood surface
    [MaxLL, i] = max(LogLikelihood( : )); % get max likelihood MaxLL and its linear index i
    [I, J] = ind2sub(szL, i); % subscripts I, J from linear index i
    
    %compute the profile likelihood functions for both parameters. (These are used
    %to compute confidence intervals).
    profileMethyFrac = max(LogLikelihood, [], 2);
    profileRate = max(LogLikelihood, [], 1);
    rateCI = -2 * (profileRate - LogLikelihood(I, J));
    methyFracCI = -2 * (profileMethyFrac - LogLikelihood(I, J));
    
    %CIs will be computed based on Likelihood Ratio Test, using Chi-squared
    %values
    perctile95ChiVal = 3.841; %95th of chi-squared distribution, 1 free param. 
    
    %Find the 95% confidence intervals
    rateRegion95 = rateGrid(rateCI <= perctile95ChiVal);
    methylFracRegion95 = methyFracGrid(methyFracCI <= perctile95ChiVal);
    rateCI95 = [rateRegion95(1) rateRegion95(end)];
    methyFracCI95 = [methylFracRegion95(1) methylFracRegion95(end)];
    
    perctile75ChiVal = 1.32; %75th of chi-squared distribution
    
    %Find the inner (75th%) confidence intervals
    rateRegion75 = rateGrid(rateCI <= perctile75ChiVal);
    methyFracCI75 = methyFracGrid(methyFracCI <= perctile75ChiVal);
    rateCI75 = [rateRegion75(1) rateRegion75(end)];
    methyFracCI75 = [methyFracCI75(1) methyFracCI75(end)];
    
    %check to see whether the inner confidence interval hits the boundary
    CheckCIsL = [minRate, maxRate] - rateCI75;
    keepRateInd = find(abs(CheckCIsL) > realmin);
    
    if numel(keepRateInd) == 0 %this indicates the inner CI overlaps entire region-->k is unidentifiable. Set k == 0
        %this was tested--the sites with only unmethylated reads are found
        %this way, and this method gives identical results to simply
        %searching for reads with only unmethylated sites
        inferedRates(ii, : ) = [0, 0, 0];
        inferredMethyFrac(ii, : ) = [0, 0, 0];
    elseif numel(keepRateInd) == 1 %this indicates CI value hits one boundary-->indicated upper/lower bound on k
        %when the CI hits the boundary, use the inner CI value as the rate
        %estimate
        inferedRates(ii, : ) = [rateCI75(keepRateInd), rateCI95(1), rateCI95(2)];
        rateInd = find(rateGrid == rateCI75(keepRateInd));
        [maxVal, maxLLind] = max(LogLikelihood( : , rateInd));
        inferredMethyFrac(ii, : ) = [methyFracGrid(maxLLind), methyFracCI95(1), methyFracCI95(2)];
    else
        inferedRates(ii, : ) = [rateGrid(J), rateCI95(1), rateCI95(2)];
        inferredMethyFrac(ii, : ) = [methyFracGrid(I), methyFracCI95(1), methyFracCI95(2)];
    end
end
save(inferredRateSavingPath, 'fittedSites', 'inferedRates', 'inferredMethyFrac')
toc
end

