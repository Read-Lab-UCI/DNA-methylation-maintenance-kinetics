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

