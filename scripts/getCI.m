function [RawRate, RawFrac, CIRate, CIFrac, MaxLL] = getCI(rateGrid, methyFracGrid, logLikelihood, prcVal)
% getCI - get confidence intervals of rate of remethylation and logLikelihood
% 
% Inputs:
%    rateGrid - the rate grid to calculate likelihood, 1 x n1 double vector
%    methyFracGrid - the methylation fraction grid to calculate likelihood, 1 x n2 double vector
%	 logLikelihood -  the likelihood surface
%	 prcVal - the percentile value of chi-squared distribution with 1 free param
% Outputs:
%    RawRate - the raw max likelihood parameters of rate
%    RawFrac - the raw max likelihood parameters of methylation fraction
%	 CIRate - confidence intervals of rate of remethylation and logLikelihood
%	 CIFrac - confidence intervals of methylation fraction
%	 MaxLL - the max value of logLikelihood

    %find the parameter values that maximize the logLikelihood
    szL=size(logLikelihood);
    [MaxLL,i]=max(logLikelihood(:));
    [I,J]=ind2sub(szL,i); %I,J are the indices of ML parameters
    
    %these are the "raw" ML parameters, but need to deal with edge cases (below)
    RawRate=rateGrid(J);
    RawFrac=methyFracGrid(I);
    
    %compute the profile likelihood function.
    ProfileRate=max(logLikelihood,[],1);
    ProfileFrac=max(logLikelihood,[],2);

    %find indices of k values with LL (logLikelihood) values within CI
    InnerInds=find(-2*(ProfileRate-logLikelihood(I,J)) <= prcVal);
    
    %Find the confidence interval
    RateRegion=rateGrid(InnerInds);
    FracRegion=methyFracGrid(-2*(ProfileFrac-logLikelihood(I,J)) <= prcVal);
    
    %write an array: [CI LB, CI UB, LL LB, LL UB]
    CIRate=[RateRegion(1) RateRegion(end) ProfileRate(InnerInds(1)) ProfileRate(InnerInds(end))];
    CIFrac=[FracRegion(1) FracRegion(end)];
end

