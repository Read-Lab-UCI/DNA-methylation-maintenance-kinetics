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

