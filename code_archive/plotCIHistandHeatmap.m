function plotCIHistandHeatmap(RateCIs, FracCIs, nrow, ncol, colID, nbins, dataName)
%plotCIHistandHeatmap - Plot the histogram of the confidence interval (CI) of methylation rates and fraction and the 2d heatmap rate CIs vs Frac CIs
%
% Inputs:
%    RateCIs - The rates to plot
%    FracCIs - The methylation fractions to plot
%	 nrow - The tot number of rows
%	 ncol - The tot number of columns
%	 colID - The column index for this data to plot
%	 nbins - The number of bins for histogram
%	 dataName - The name of input data
	
	subplot(nrow, ncol, colID);
	histogram(RateCIs, nbins, 'Normalization', 'probability');
	ylabel('Probability');
	% xlim([0, 10]);
	ylim([0, 0.7]);
	title(strcat('(RCI/Rate): ',dataName));

	subplot(nrow, ncol, 1 * ncol + colID);
	histogram(FracCIs, nbins, 'Normalization', 'probability');
	ylabel('Probability');
	xlim([0, 1]);
	ylim([0, 0.25]);
	title(strcat('CI of Frac: ',dataName));

	subplot(nrow, ncol, 2 * ncol +colID);

	rateCIGrid = 10.^(-1: .1 : 3); % grid of the CI for rate-values
    methyFracCIGrid = 0: .1 : 1.0; %grid of the CI of f-values

	heatmapDataArray=zeros(numel(rateCIGrid), numel(methyFracCIGrid));
    for ii = 1 : numel(RateCIs)
        indsa = find(RateCIs(ii) <= rateCIGrid);
        indsb = find(FracCIs(ii) <= methyFracCIGrid);
        heatmapDataArray(indsa(1),indsb(1)) = heatmapDataArray(indsa(1),indsb(1)) + 1;
    end
    heatmapDataArray

    heatmapDataArray=heatmapDataArray'/numel(RateCIs);

	CC = corrcoef(RateCIs, FracCIs) % Correlation coefficient between methylation rates and methylation fraction
	logRateCIGrid = log10(rateCIGrid);
	
	imagesc(logRateCIGrid, methyFracCIGrid, -log(heatmapDataArray));
 	colormap(flipud(parula));
 	h = colorbar;
 	Vals = fliplr([0.0001, 0.001, 0.01, 0.1]);
 	Ticks = -log(Vals);
 	set(h,'YDir', 'reverse');
 	set(h,'Ticks', Ticks);
 	set(h,'TickLabels', cellstr(num2str(Vals(:)))');
 	set(gca, 'YDir', 'normal');
 	axis square
 	xlabel('log10(RCI/Rate)');
 	if (colID == 1)
		ylabel('CI of Frac');
 		ylabel(h, 'Probability');
 		title(strcat('#: ', int2str(length(RateCIs)), ', ', dataName, ', Corr: ', num2str(round(CC(2, 1), 2))));
 	else
 		title(strcat(dataName, ', Corr: ', num2str(round(CC(2, 1), 2))));
 	end

end