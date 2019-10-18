function plotConfidenceIntervalHistogram(inferredRatePath, figSavingDir)
%plotConfidenceIntervalHistogram - Plot the histogram of the confidence interval (CI) of methylation rates and fraction
%
% Inputs:
%    inferredRatePath - Path of the inferred kinetic rates and methylation fraction
%    figSavingDir - Directory to save figures
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

	rateCIs = inferedRates(:, 3) - inferedRates(:, 2) ; % The CIs for rates
	fracCIs = inferredMethyFrac(:, 3) - inferredMethyFrac(:, 2) ; % The CIs for methylaiton fractions

	%CI histogram of fitted methylation rates
	slowRateThreshold = prctile(Rates, 33.3); % the slow kinetic rate(k): if k < 33.3 percentile
	fastRateThreshold = prctile(Rates, 66.6); % the fast kinetic rate(k): if k > 66.6 percentile

	slowIndexs = find(Rates > 0.0 & Rates <= slowRateThreshold); % the indexs of slow rates
	slowRateCIs = rateCIs(slowIndexs); % the rates CIs of slow rates
	slowFracionCIs = fracCIs(slowIndexs) ; % the methylation fractions CIs of slow rates
	
	midiumIndexs = find(Rates > slowRateThreshold & Rates < fastRateThreshold );
	midiumRateCIs = rateCIs(midiumIndexs); 
	midiumFracionCIs = fracCIs(midiumIndexs) ;

	fastIndexs = find(Rates >= fastRateThreshold);
	fastRateCIs = rateCIs(fastIndexs);
	fastFracionCIs = fracCIs(fastIndexs);

	figure(1)
	
	fontSize = 20;  %font size of figure
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 6])
	set(gca, 'FontSize', fontSize);
	nbins = 20;

	maxP1 = 0.7; % max probability of histogram
	subplot(2, 4, 1);
	histogram(log10(rateCIs(find(Rates > 0.0))), nbins, 'Normalization', 'probability');
	xlim([-2, 1.1]);
	ylabel('Probability');
	ylim([0, maxP1]);
	title('CI of All Rates');
	
	subplot(2, 4, 2);
	histogram(log10(slowRateCIs), nbins, 'Normalization', 'probability');
	xlim([-2, 1.1]);
	ylim([0, maxP1]);
	xlabel('CI of Kinetic Rates');
	title('CI of Slow Rates');

	subplot(2, 4, 3);
	histogram(log10(midiumRateCIs), nbins, 'Normalization', 'probability');
	xlim([-2, 1.1]);
	ylim([0, maxP1]);
	title('CI of Midium Rates');
	
	subplot(2, 4, 4);
	histogram(log10(fastRateCIs), nbins, 'Normalization', 'probability');
	xlim([-2, 1.1]);
	ylim([0, maxP1]);
	title('CI of Fast Rates');

	maxP2 = 0.25;
	subplot(2, 4, 5);
	histogram(fracCIs(find(Rates > 0.0)), nbins, 'Normalization', 'probability');
	ylabel('Probability');
	ylim([0, maxP2]);
	title('Frac CI');
	
	subplot(2, 4, 6);
	histogram(slowFracionCIs, nbins, 'Normalization', 'probability');
	ylim([0, maxP2]);
	xlabel('CI of Methylation Fraction');
	title('of Slow Rates');

	subplot(2, 4, 7);
	histogram(midiumFracionCIs, nbins, 'Normalization', 'probability');
	ylim([0, maxP2]);
	title('of Midium Rates');

	subplot(2, 4, 8);
	histogram(fastFracionCIs, nbins, 'Normalization', 'probability');
	ylim([0, maxP2]);
	title('of Fast Rates');
	
	fig1Path = strcat(figSavingDir, "ConfidenceIntervalHistogram.eps");
	print(fig1Path, '-depsc', '-r300')
	close all

end