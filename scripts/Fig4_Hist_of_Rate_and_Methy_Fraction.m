function fig4HistofRateandMethyFraction(inferredRatePath, figSavingDir)
%fig4HistofRateandMethyFraction - Plot the histogram of kinetic rates and steady state methylation fraction, separately (Fig 4 in paper)
%
% Inputs:
%    inferredRatePath - Path of the inferred kinetic rates and methylation fraction
%    figSavingDir - Directory to save figures
%
% Example: 
%    fig4HistofRateandMethyFraction('../data/InferedRates/kineticRateChr1.mat', '../figures')

	% if figSavingDir not exist, make a directory
	if ~exist(figSavingDir)
    	mkdir(figSavingDir);
	end

	load(inferredRatePath, 'fittedSites', 'inferedRates', 'inferredMethyFrac');

	Rates = inferedRates(:, 1);
	Fracs = inferredMethyFrac(:, 1);

	%Plot param configuration
	
	fontSize = 20;%font size of figure

	%histogram of fitted methylation rates
	figure(1)
	h1 = histogram(log10(Rates), 'Normalization', 'probability', 'BinWidth', 0.15)
	set(gca, 'FontSize', fontSize)
	xlim([-2.2 1.5])
	xlabel('log10 Rate Constant k (hr^{-1})')
	ylabel('Probability')
	title('Distribution of Fitted Methylation Rates');
	fig1Path = strcat(figSavingDir, "Fig4a_rate_histogram.eps");
	print(fig1Path, '-depsc')

	%histogram of fitted steady state methylation fractions
	figure(2)
	h2 = histogram(Fracs, 'Normalization', 'probability', 'BinWidth', 0.05)
	set(gca, 'FontSize', fontSize)
	title('Distribution of Fitted Fraction Methylation')
	xlabel('Fraction of Cells with Methylation/per Site (f)')
	ylabel('Probability')
	fig2Path = strcat(figSavingDir, "Fig4b_methyfraction_histogram.eps");
	print(fig2Path, '-depsc')
	close all
end

