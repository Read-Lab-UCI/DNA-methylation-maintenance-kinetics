function plotHeatmapforInferredRateandFraction(inferredRatePath, figSavingDir)
%plotHeatmapforInferredRateandFraction - Plot the heatmap for kinetic rates and steady state methylation fraction
%
% Inputs:
%    inferredRatePath - Path of the inferred kinetic rates and methylation fraction
%    figSavingDir - Directory to save figures
%
% Example: 
%    plotHeatmapforInferredRateandFraction('../data/InferedRates/kineticRateChr1.mat', '../figures')

	% if figSavingDir not exist, make a directory
	if ~exist(figSavingDir)
    	mkdir(figSavingDir);
	end

	load(inferredRatePath, 'fittedSites', 'inferedRates', 'inferredMethyFrac');

	Rates = inferedRates(:, 1);
	Fracs = inferredMethyFrac(:, 1);
	
	logMinRate = -2; %log10 minimum allowed fit rate
    logMaxRate = 1;  %log10 maximume allowed fit rate
    rateGrid = 10.^(logMinRate : .1 : logMaxRate); % grid of k-values
    methyFracGrid = 0 : 0.1 : 1; %grid of f-values at which LogLikelihood will be computed

    heatmapDataArray=zeros(numel(rateGrid), numel(methyFracGrid));
    for ii = 1 : numel(fittedSites)
        indsa = find(Rates(ii) <= rateGrid);
        indsb = find(Fracs(ii) <= methyFracGrid);
        heatmapDataArray(indsa(1),indsb(1)) = heatmapDataArray(indsa(1),indsb(1)) + 1;
    end


	heatmapDataArray=heatmapDataArray/numel(fittedSites);
	
	CC = corrcoef(Rates, Fracs) % Correlation coefficient between methylation rates and methylation fraction
	logRateGrid = log10(rateGrid);
	
	f = imagesc(logRateGrid, methyFracGrid, -log(heatmapDataArray));
	colormap(flipud(parula));
	h = colorbar;
	ylabel(h, 'Probability');
	Vals = fliplr([0.0001, 0.001, 0.01, 0.1]);
	Ticks = -log(Vals);
	set(h,'YDir', 'reverse');
	set(h,'Ticks', Ticks);
	set(h,'TickLabels', cellstr(num2str(Vals(:)))');
	set(gca, 'YDir', 'normal');
	set(gca, 'FontSize', 20);
	axis square
	xlabel('log10 Rate Constant k (hr^{-1})');
	ylabel('Methy Fraction');
	title(strcat('Rate vs Fraction, Correlation: ', num2str(round(CC(2,1), 2))));
	figPath = strcat(figSavingDir, "Heatmap_For_Inferred_Rate_and_Fraction.eps");
	print(figPath, '-depsc');
	close all
end

