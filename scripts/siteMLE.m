function [inferedRateNCIs, inferredMethyFracNCIs, LogLikelihood] = siteMLE(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii)
% calcLikelihoodSurface - calculate the likelihood surface
% 
% Inputs:
%    rateGrid - the rate grid to calculate likelihood, 1 x n1 double vector
%    methyFracGrid - the methylation fraction grid to calculate likelihood, 1 x n2 double vector
%    timePoints - the timepoints at which the reads are experimentally measured, 1 x 4 double vector
%    methylatedReadTimeCourse - the time course counts of methylated reads for the given timepoints, 1 x 4 double vector
%    unmethylatedReadTimeCourse - the time course counts of unmethylated reads for the given timepoints, 1 x 4 double vector
% Outputs:
%    inferedRateNCIs - inferred kinetic rates and its confidence interval
%    inferredMethyFracNCIs - inferred methylation fraction and its confidence interval
%
% Example: 
%    [inferedRateNCIs, inferredMethyFracNCIs] = siteMLE(10.^(-2: .01: 1), 0: 0.01: 1, [0.5, 1.5, 4.5, 16.5], [7, 0, 2, 3], [1, 0, 5, 6])
    
    LogLikelihood = calcLikelihoodSurface(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii);
    
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
    maxRate = max(rateGrid);
    minRate = min(rateGrid);
    CheckCIsL = [minRate, maxRate] - rateCI75;
    keepRateInd = find(abs(CheckCIsL) > realmin);
    
    if numel(keepRateInd) == 0 %this indicates the inner CI overlaps entire region-->k is unidentifiable. Set k == 0
        %this was tested--the sites with only unmethylated reads are found
        %this way, and this method gives identical results to simply
        %searching for reads with only unmethylated sites
        inferedRateNCIs = [0, 0, 0];
        inferredMethyFracNCIs = [0, 0, 0];
    elseif numel(keepRateInd) == 1 %this indicates CI value hits one boundary-->indicated upper/lower bound on k
        %when the CI hits the boundary, use the inner CI value as the rate
        %estimation
        
        %if k is on boundary (lower or upper)
        %copy above code (lines 104-110) to introduce X% where X is desired (low) threshold
        %calculate the X% threshold. 
        %increase the k-domain boundary (or decrease it if k is on slowest
        %lower boundary. 

        %Recalculate the log-likelihood surface and recalculate
        %X threshold and 95 threshold. If abs(CI_X_new-CI_X_old)< epsilon (epsilon = 0.1) then stop
        %use abs(CI_X-CI_95) as the confidence interval. 
        %report: for right boundary: inferredRates(ii,:)=[rateCI75(keepRateInd),rateCI95(1),
        %rateCIX(1)]; (or switch order for left boundary). 
    
        inferedRateNCIs = [rateCI75(keepRateInd), rateCI95(1), rateCI95(2)];
        rateInd = find(rateGrid == rateCI75(keepRateInd));
        [maxVal, maxLLind] = max(LogLikelihood( : , rateInd));
        inferredMethyFracNCIs = [methyFracGrid(maxLLind), methyFracCI95(1), methyFracCI95(2)];
    else
        inferedRateNCIs = [rateGrid(J), rateCI95(1), rateCI95(2)];
        inferredMethyFracNCIs = [methyFracGrid(I), methyFracCI95(1), methyFracCI95(2)];
    end
end

