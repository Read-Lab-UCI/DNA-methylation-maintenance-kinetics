function [] = MLEInference(AllDatDir, RatesDir, StartChr, EndChr)
    clc;
    close all;
    if ~exist(RatesDir, 'dir')
       mkdir(RatesDir)
    end

    EnoughReads0=10; %number of reads required at t=0
    EnoughReadsLater=5; %number of reads required for later timepoints
    
    for chromosome = StartChr : EndChr
        inputReadDataPath = strcat(AllDatDir,'/AllDat_chr',int2str(chromosome),'.mat'); % file path of read data
        inferredRateSavingPath = strcat(RatesDir, '/Rates_chr',int2str(chromosome),'.mat');
        [fittedSites, inferedRates, inferredMethyFrac, AllDat] = kineticRatesMLE(inputReadDataPath,inferredRateSavingPath, EnoughReads0, EnoughReadsLater);
    end
end

function [fittedSites, inferedRates, inferredMethyFrac, AllDat] = kineticRatesMLE(inputReadDataPath,inferredRateSavingPath, EnoughReads0, EnoughReadsLater)
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
disp(sprintf('Inferring for %s now ...', inputReadDataPath));
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

Reads=sum(AllDat(:,:,1:2),3);
NumReads0=Reads(:,1);
NumReadsLater=sum(Reads(:,2:end),2);
KeepSites=find(NumReadsLater>=EnoughReadsLater & NumReads0>=EnoughReads0);
numSites=numel(KeepSites);

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
fittedSites = sites(KeepSites)';
AllDat = AllDat(KeepSites, : , :);
% siteidx = find(fittedSites == 91267)
% fittedSites = fittedSites(siteidx);
tic

%loop over the sites for rate and steady state methylation level MLE
parfor ii = 1 : numSites  %parfor
    %get the read data for site ii
    methylatedReadTimeCourseForSiteii = AllDat(ii, : , 1);
    unmethylatedReadTimeCourseForSiteii = AllDat(ii, : , 2);
    [inferedRates(ii, : ), inferredMethyFrac(ii, : ), LogLikelihood] = siteMLE(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii);
%     inferedRates(ii, : )
%     inferredMethyFrac(ii, : )
end
save(inferredRateSavingPath, 'fittedSites', 'inferedRates', 'inferredMethyFrac')
toc
end

function [inferedRateAndCIs, inferredMethyFracAndCIs, LogLikelihood] = siteMLE(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii)
% calcLikelihoodSurface - calculate the likelihood surface
% 
% Inputs:
%    rateGrid - the rate grid to calculate likelihood, 1 x n1 double vector
%    methyFracGrid - the methylation fraction grid to calculate likelihood, 1 x n2 double vector
%    timePoints - the timepoints at which the reads are experimentally measured, 1 x 4 double vector
%    methylatedReadTimeCourse - the time course counts of methylated reads for the given timepoints, 1 x 4 double vector
%    unmethylatedReadTimeCourse - the time course counts of unmethylated reads for the given timepoints, 1 x 4 double vector
% Outputs:
%    inferedRateAndCIs - inferred kinetic rates and its confidence interval
%    inferredMethyFracAndCIs - inferred methylation fraction and its confidence interval
%	 LogLikelihood - the likelihood surface
% Example: 
%    [inferedRateAndCIs, inferredMethyFracAndCIs] = siteMLE(10.^(-2: .01: 1), 0: 0.01: 1, [0.5, 1.5, 4.5, 16.5], [7, 0, 2, 3], [1, 0, 5, 6])
    
    LogLikelihood = calcLikelihoodSurface(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii);
    
    %CIs will be computed based on Likelihood Ratio Test, using Chi-squared values
    PrcVal95=3.841; %95th of chi-squared distribution, 1 free param.
    PrcVal75=1.32; %75th
    PrcValEdge=0.000157; %1th

    maxRate = max(rateGrid);
	minRate = min(rateGrid);
	
	logMinRate = log10(maxRate);
	logMaxRate = log10(maxRate);

    %Find the edge confidence intervals
    [RawRate,RawFrac,CIRate95,CIfrac95,MaxLL]= getCI(rateGrid, methyFracGrid, LogLikelihood, PrcVal95);

    CheckCIsL=[minRate, maxRate]-CIRate95(1:2);
    keepIndRate=find(abs(CheckCIsL)>realmin);

    if numel(keepIndRate) == 0 %this indicates the CI overlaps entire region-->k is unidentifiable. Set k==0
        %this was tested--the sites with only unmethylated reads or very few methylated reads are found
        inferedRateAndCIs = [0, 0, 100, -1];
        inferredMethyFracAndCIs = [0, 0, 1];
    elseif abs(maxRate - CIRate95(2)) <= realmin %if upper CI 95 lies on boundary
        %loop to keep increasing upper k boundary to see where LL stops decreasing
        LLCIRate95_Old=CIRate95(4);
        Diff_Old = 1;
        while Diff_Old > PrcValEdge
            logMaxRate = logMaxRate + 0.05; %increase max k
            rateGrid = [rateGrid, 10^logMaxRate]; %expanded grid of k values
            maxRate = max(rateGrid);
            LogLikelihood=calcLikelihoodSurface(rateGrid,methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii);
            [RawRate,RawFrac,CIRate95]= getCI(rateGrid, methyFracGrid, LogLikelihood, PrcVal95);
            Diff_Old=abs(CIRate95(4)-LLCIRate95_Old);
            LLCIRate95_Old=CIRate95(4);
        end

        %recompute CIs, in case they have changed with expanded parameter grid
        [RawRate,RawFrac,CIRate95,CIfrac95]= getCI(rateGrid, methyFracGrid, LogLikelihood, PrcVal95);
        if abs(maxRate - RawRate) <= realmin %if k lies on boundary (only lower k limit may be determined)
            %set the k estimates by: Estimated Value: 75th LB. Note
            %that the 95th UB is only an estimate--a lower bound on the
            %position of 95th upper CI (since the data do not constrain
            %the upper limit)
            [RawRate,RawFrac,CIRate75]= getCI(rateGrid, methyFracGrid, LogLikelihood,PrcVal75);
            
            RateEst=CIRate75(1);
            %now find the corresponding max f value
            [aa, bb]=max(LogLikelihood(:,find(rateGrid==RateEst)));
            inferredMethyFracAndCIs=[methyFracGrid(bb),CIfrac95(1),CIfrac95(2)];
            inferedRateAndCIs=[RateEst,CIRate95(1),CIRate95(2),1];
        else
            inferedRateAndCIs=[RawRate,CIRate95(1),CIRate95(2),2];
            inferredMethyFracAndCIs=[RawFrac,CIfrac95(1),CIfrac95(2)];
        end
    elseif abs(minRate-CIRate95(2))<=realmin %if lower 95 CI lies on boundary loop to keep decreasing lower k boundary to see where LL stops decreasing
        LLCIRate95_Old=CIRate95(3);
        Diff_Old=1;
        while Diff_Old>PrcValEdge
            logMinRate = logMinRate-0.05; %decrease min k
            rateGrid=[10^logMinRate,rateGrid]; %expanded grid of k values
            minRate=min(rateGrid);
            LogLikelihood=calcLikelihoodSurface(rateGrid,methyFracGrid, timePoints, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii);
            [RawRate,RawFrac,CIRate95]= GetCILam_Func(rateGrid,methyFracGrid,LogLikelihood,PrcValEdge);
            Diff_Old=abs(CIRate95(3)-LLCIRate95_Old);
            LLCIRate95_Old=CIRate95(3);
        end
        %recompute CIs, in case they have changed with expanded parameter grid
        [RawRate,RawFrac,CIRate95]= GetCILam_Func(rateGrid,methyFracGrid,LogLikelihood,PrcVal95);
        [RawRate,RawFrac,CIfrac95]= GetCIFrac_Func(rateGrid,methyFracGrid,LogLikelihood,PrcVal95);
        if abs(minRate-RawRate)<=realmin %if k lies on boundary
            %set the k estimates by: Estimated Value: 75th UB. CI LB: 1th
            %UB, CI UB: 95th upper
            [RawRate,RawFrac,CIRate75]= GetCILam_Func(rateGrid,methyFracGrid,LogLikelihood,PrcVal75);
            RateEst=CIRate75(2);
            %now find the corresponding max f value
            [aa,bb]=max(LogLikelihood(:,find(rateGrid==CIRate75(2))));
            inferredMethyFracAndCIs=[methyFracGrid(bb),CIfrac95(1),CIfrac95(2)];
            %LB=max([0,CIRate95(1)]);
            inferedRateAndCIs=[RateEst,CIRate95(1),CIRate95(2),1];
        else
            %LB=max([0,CIRate95(1)]);
            inferedRateAndCIs=[RawRate,CIRate95(1),CIRate95(2),2];
            inferredMethyFracAndCIs=[RawFrac,CIfrac95(1),CIfrac95(2)];
        end
    else %if no edge case is relevant
        inferedRateAndCIs=[RawRate,CIRate95(1),CIRate95(2),0];
        inferredMethyFracAndCIs=[RawFrac,CIfrac95(1),CIfrac95(2)];
    end
end


function [LogLikelihood] = calcLikelihoodSurface(rateGrid, methyFracGrid, timePoints, methylatedReadTimeCourse, unmethylatedReadTimeCourse)
% calcLikelihoodSurface - calculate the likelihood surface
% 
% Inputs:
%    rateGrid - the rate grid to calculate likelihood, 1 x n1 double vector
%    methyFracGrid - the methylation fraction grid to calculate likelihood, 1 x n2 double vector
%    timePoints - the timepoints at which the reads are experimentally measured, 1 x 4 double vector
%    methylatedReadTimeCourse - the time course counts of methylated reads for the given timepoints, 1 x 4 double vector
%    unmethylatedReadTimeCourse - the time course counts of unmethylated reads for the given timepoints, 1 x 4 double vector
% Outputs:
%    LogLikelihood - The logLikelihood surface for given rate and methylation fraction grid, n1 x n2 double matrix
%
% Example: 
%    calcLikelihoodSurface(10.^(-2: .01: 1), 0: 0.01: 1, [0.5, 1.5, 4.5, 16.5], [7, 0, 2, 3], [1, 0, 5, 6])

    pMethyArray = zeros(numel(methyFracGrid), numel(rateGrid), numel(timePoints));
    pUnmethyArray = zeros(numel(methyFracGrid), numel(rateGrid), numel(timePoints));

    %loop over the timepoints to compute the LogLikelihood surface
    for tind = 1 : numel(timePoints)
        timeI = timePoints(tind);
        expoentialTerms = 1 - exp(-rateGrid .* timeI);
        expoentialTerms(expoentialTerms >= 1) = 1 - eps;
        repExpoentialTerms = repmat(expoentialTerms, numel(methyFracGrid), 1);
        repMethyFrac = repmat(methyFracGrid', 1, numel(rateGrid));
        pMethyArray( : , : , tind) = methylatedReadTimeCourse(tind) .* log(repMethyFrac .* repExpoentialTerms);
        pUnmethyArray( : , : , tind) = unmethylatedReadTimeCourse(tind) .* log(1 - repMethyFrac .* repExpoentialTerms);
    end
    %this is the LogLikelihood surface as a function of the parameters
    LogLikelihood = sum(pMethyArray, 3, 'omitnan') + sum(pUnmethyArray, 3, 'omitnan');
end


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

