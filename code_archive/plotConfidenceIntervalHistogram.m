function plotConfidenceIntervalHistogram(inferredRatePath, figSavingDir, maxRateCI, maxFracCI)
%plotConfidenceIntervalHistogram - Plot the histogram of the confidence interval (CI) of methylation rates and fraction
%
% Inputs:
%    inferredRatePath - Path of the inferred kinetic rates and methylation fraction
%    figSavingDir - Directory to save figures
%	 maxRateCI - Max CI value allowed for re-methylaiton rates 
%	 maxFracCI - Max CI value allowed for methylation fraction
%
% Example: 
%    plotConfidenceIntervalHistogram('../data/InferedRates/kineticRateChr1.mat', '../figures')

	% if figSavingDir not exist, make a directory
	if ~exist(figSavingDir)
    	mkdir(figSavingDir);
	end

	load(inferredRatePath, 'fittedSites', 'inferedRates', 'inferredMethyFrac');
	Rates = inferedRates(:, 1);
	Fracs = inferredMethyFrac(:, 1);

	rateCIs = (inferedRates(:, 3) - inferedRates(:, 2)) ./ Rates; % The CIs for rates
	fracCIs = inferredMethyFrac(:, 3) - inferredMethyFrac(:, 2); % The CIs for methylaiton fractions
	allIndexs = find((Rates > 0.0) & (rateCIs < maxRateCI)   & (fracCIs < maxFracCI));

	allRates = Rates(allIndexs);
	allFracs = Fracs(allIndexs);
	allRateCIs = rateCIs(allIndexs);
	allFracCIs = fracCIs(allIndexs);

	%CI histogram of fitted methylation rates
	slowRateThreshold = prctile(Rates, 33.3); % the slow kinetic rate(k): if k < 33.3 percentile
	fastRateThreshold = prctile(Rates, 66.6); % the fast kinetic rate(k): if k > 66.6 percentile

	slowIndexs = find((Rates > 0.0) & (Rates <= slowRateThreshold  )  & (rateCIs < maxRateCI)   & (fracCIs < maxFracCI)); % the indexs of slow rates
	slowRateCIs = rateCIs(slowIndexs); % the rates CIs of slow rates
	slowRates = Rates(slowIndexs);
	slowFracs = Fracs(slowIndexs);
	slowFracionCIs = fracCIs(slowIndexs) ; % the methylation fractions CIs of slow rates
	
	midiumIndexs = find((Rates > slowRateThreshold) & (Rates < fastRateThreshold )  & (rateCIs < maxRateCI)   & (fracCIs < maxFracCI));
	midiumRates = Rates(midiumIndexs);
	midiumFracs = Fracs(midiumIndexs);
	midiumRateCIs = rateCIs(midiumIndexs); 
	midiumFracionCIs = fracCIs(midiumIndexs) ;

	fastIndexs = find((Rates >= fastRateThreshold)  & (rateCIs < maxRateCI)   & (fracCIs < maxFracCI));
	fastRates = Rates(fastIndexs);
	fastFracs = Fracs(fastIndexs);
	fastRateCIs = rateCIs(fastIndexs);
	fastFracionCIs = fracCIs(fastIndexs);

	figure(1)
	
	fontSize = 20;  %font size of figure
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
	set(gca, 'FontSize', fontSize);
	nbins = 20;
	nrow = 5;
	ncol = 4;
	plotCIHistandHeatmap(allRateCIs, allFracCIs, allRates, allFracs, nrow, ncol, 1, nbins, 'All Rates');
	plotCIHistandHeatmap(slowRateCIs, slowFracionCIs, slowRates, slowFracs, nrow, ncol, 2, nbins, 'Slow Rates');
	plotCIHistandHeatmap(midiumRateCIs, midiumFracionCIs, midiumRates, midiumFracs, nrow, ncol, 3, nbins, 'Midium Rates');
	plotCIHistandHeatmap(fastRateCIs, fastFracionCIs, fastRates, fastFracs, nrow, ncol, 4, nbins, 'Fast Rates');

	fig1Path = strcat(figSavingDir, "ConfidenceIntervalHistogram.png");
	print(fig1Path, '-dpng', '-r300')
	close all

end