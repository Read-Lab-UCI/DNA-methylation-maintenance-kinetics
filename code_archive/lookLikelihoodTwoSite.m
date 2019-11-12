function [] = lookLikelihoodTwoSite()
    close all
    clear
    tShift = 0.5;
    Times = [0, 1, 4, 16] + tShift;

    fracgrid = 0 : 0.01 : 1;    %grid of f-values at which LogLikelihood will be computed
    lamgrid = 10.^(-2 : .01 : 1);   %grid of k-values

    inputReadDataPath = '../data/ReadData/AllDat_Chr1WT.mat';
    load(inputReadDataPath,'AllDat','sites');
    Reads = sum(AllDat(:, :, 1:2), 3);
    NumReads = sum(Reads, 2);
    
    inferredRatePath = '../data/InferedRates/kineticRateChr1.mat';
    load(inferredRatePath,'inferedRates', 'inferredMethyFrac');
    Rates = inferedRates(:, 1);

    HighKLowReads = find(Rates  > 5 & NumReads < 10);
    HighKHighReads = find(Rates > 5 & NumReads > 20);
    
    MethReads=AllDat(:,:,1);
    NumMethReads=sum(MethReads,2);
    NumTimPts=sum(Reads > 0, 2);
    Tots = sum(AllDat,3);
    PercMethSites = AllDat(:,:,1)./Tots;

    siteIndex = 20;
    Listinds= [HighKLowReads(siteIndex), HighKHighReads(siteIndex)]
    FS = 10;
    figure(1);
    ha=tight_subplot(2,3,[0.08 0.07],[0.09 0.01],[.09 0.01]);
        
    cd ../scripts
    for loopind=1:numel(Listinds)
        ind=Listinds(loopind);
        
        methylatedReadTimeCourseForSiteii = AllDat(ind,:,1);
        unmethylatedReadTimeCourseForSiteii = AllDat(ind,:,2);
        
        [MLELam, MLEFrac, LogLikelihood] = siteMLE(lamgrid, fracgrid, Times, methylatedReadTimeCourseForSiteii, unmethylatedReadTimeCourseForSiteii);

        ProfileLam = max(LogLikelihood, [], 1);
        ProfileFrac = max(LogLikelihood, [], 2);

        NormLam = exp(ProfileLam-max(ProfileLam));
        NormFrac = exp(ProfileFrac-max(ProfileFrac));

        rowind=(loopind-1)*3;
        axes(ha(rowind+1))
        plot(lamgrid,NormLam,'-k','LineWidth',2)
        %semilogx(lamgrid,NormLam,'-k','LineWidth',2)
        hold on
        yvals=[0 1.2*max(NormLam)];
        xvals=repmat(MLELam(2:3),numel(yvals),1);
        plot(xvals(:,1),yvals,'-.r','LineWidth',1.5)

        if xvals(1, 2) < lamgrid(end)
            plot(xvals(:,2),yvals,'-.r','LineWidth',1.5)
        end

        xlim([min(lamgrid)/1.2 max(lamgrid)*1.2])
        set(gca,'XTick',[2 5 7 10])
        ylim([0 1.05])
        
        if loopind==numel(Listinds)
             xlabel('Methylation Rate (k)')
        end
        if loopind==1
             legend('Likelihood','95% C.I.','Location','West')
        end

        ylabel('Relative Likelihood')
        set(gca,'FontSize',FS)
        
        axes(ha(rowind+2))
        plot(fracgrid, NormFrac,'-k','LineWidth',2)
        hold on
        yvals=[0 1.2 * max(NormFrac)];
        xvals=repmat(MLEFrac(2:3),numel(yvals),1);
        plot(xvals(:,1),yvals,'-.r','LineWidth',1.5)
        if xvals(1,2)<fracgrid(end)
            plot(xvals(:,2),yvals,'-.r','LineWidth',1.5)
        end
        if loopind==numel(Listinds)
             xlabel('Methylation Fraction (f)')
        end
        xlim([-0.05 1])
        ylim([0 1.05])
        
        set(gca,'FontSize',FS)
        
        ModelTimes=[0:.1:17];
        Model=MLEFrac(1)*(1-exp(-MLELam(1)*ModelTimes));
        
        NSamples=10000;
        Samples=cell(1,4);
        for tind=1:numel(Times)
            for ns=1:NSamples
                %Meth0=find(rand(Reads(ind,tind))<=MLEFrac(1));
                PMethAlready=MLEFrac(1)*(1-exp(-MLELam(1)*Times(tind)));
                Samples{tind}=rand(Reads(ind, tind),1)<=PMethAlready;
                PercMethModel(ns,tind)=mean(Samples{tind});
                %ci(:,tind)=bootci(500,@mean,PercMethModel(:,tind));
            end
        end
        ci(1,:)=prctile(PercMethModel,2.5);
        ci(2,:)=prctile(PercMethModel,97.5);
        
        axes(ha(rowind+3))
        plot(Times,PercMethSites(ind,:),'o','Color',[0,0,1],'MarkerFaceColor','b','MarkerSize',5)
        hold on
        plot(ModelTimes, Model,'-k','LineWidth',2)
        ylim([0 1.1])
        xlim([-1 18])
       
       
        hold on
        plot(Times,ci(1,:),'-.r','LineWidth',1.5)
        plot(Times,ci(2,:),'-.r','LineWidth',1.5)
    if loopind==1
        legend('Mean Exp. Data','Mean Model','95% C.I. Model','Location','SouthEast')
    end   
        ylabel('Frac. Methylation (Nascent)')
        xlabel('Time (hrs)')
        set(gca,'FontSize',FS)
        
    end
    print -depsc2 LLExamples_TwoSite
    print -dpng LLExamples_TwoSite
end







