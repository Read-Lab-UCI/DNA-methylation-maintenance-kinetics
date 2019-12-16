function [SimDat] = Processive_FPTKMC_1DDiff(sites,fracmeth,TimePt,TrajFlag,SimFlag)
%sites and fracmeth: the nucleotide location and associated steady-state fraction methylated: these are used to intialize the time-0 state of DNA at start of simulation

%TimePt is the value of time (in hours) at which state of system is
%requested (as close as possible to requested time with Gillespie
%algorithm)

%TrajFlag: 0 to just save the requested timepoints, 1 to save the entire trajectory information
%SimFlag: 0 to "start from scratch" including expensive computation of FPT lookup tables for FPTKMC method. If these have already been computed and saved to disk, use SimFlag=1

%diffusion parameters
D=1E6*3600; %bp^2/h. Obtained diffusion coefficients for 1D sliding of transcription factors, such as LacI and p53 as reported by Mirny et al.
MaxDist=70;
latt=1; %lattice spacing in b.p.
koff=5E6;
LookupFName=['D' num2str(D) '_koff' num2str(koff) '_MD' num2str(MaxDist)];
CDFName=['CDFs' LookupFName];
FPT_tName=['FPT_ts' LookupFName];
ExitName=['Exit_Probs' LookupFName];

if SimFlag==0
    [CDFs,FPT_ts,Exit_Probs]=MakeFPTLookup(D,koff,latt,MaxDist);
    disp('table done')
    save(CDFName,'CDFs','-v7.3')
    save(FPT_tName,'FPT_ts','-v7.3')
    save(ExitName,'Exit_Probs','-v7.3')
else
    CDFs=load(CDFName,'CDFs');
    FPT_ts=load(FPT_tName,'FPT_ts');
    Exit_Probs=load(ExitName,'Exit_Probs');
end

NumCpGSites=length(sites);
NumRxns=5;  %number of reactions
NumSpec=6; % [u,h,Eh,EAh,Em,m] %species in the model
%xu is a non-CpG site that is unbound by E,
%xb is a non-CpG site that is bound by E,
%the other species are CpGs in various states of methylation/E-bound/unbound


%COTPM adapted from Pradhan
%E + h <-> Eh k1f, k1r
%EA + h <-> EAh (I added these reactions) k2f, k2r

%Processive mechanism from Svedruzic + diffusion ideas from Erban
%EQm + h -> m + Eh     %DNMT1 hop to the right (k_hop_R)
%h + EQm -> Eh + m     %DNMT1 hop to the left (k_hop_L)

%EQmi-> m koff off rates

%following calculations are based on assumption of nuclear volume=1000
%microns-cubed=1E-12 Liters, and [DNMT1]=5000 nM-500 microMolar.
V_nuc=1E-12;
h=28E6; % # CpG sites
Ratio=1/100;
h0=h/V_nuc/6.022E23*1E6;%initial [h] in micromolar
E0=h0*Ratio; %total number of DNMT1 enzymes in system
kcat=40; %/hr
Kmh=1.3;%uM
KmA=5.5;%uM
A0=100; %uM
Kia=2.85; %uM

k1f=kcat/Kmh;
k2r=100; %arbitrary value(/hr)
k1r=Kia*k1f;
k2f=(kcat+k2r)/(KmA);

E=round(E0*1E-6*6.022E23*NumCpGSites*V_nuc/h) % # of enzyme copies scaled by NumCpGSites
A=round(A0*1E-6*6.022E23*NumCpGSites*V_nuc/h); % # of AdoMet copies scaled by NumCpGSites

k1f=k1f*h/NumCpGSites/V_nuc*1E6/6.02E23; % /(copies*h) scaled by NumCpGSites
k2f=k2f*h/NumCpGSites/V_nuc*1E6/6.02E23; % /(copies*h) scaled by NumCpGSites

SimDat=zeros(NumCpGSites,2); %initialize AllDat with 4 columns, 2nd dimension is # methylated, #
%unmethylated

%setting up reaction stoichiometries for Kinetic Monte Carlo (Gillespie
%part). (The diffusive superhops are taken care of elsewhere
Stoich=zeros(NumRxns,NumSpec); %stoichiometry array
%Rxn 1; h + E -> Eh, k1f
Stoich(1,2)=-1;
Stoich(1,3)=1;
%Rxn 2; Eh -> h + E, k1r
Stoich(2,2)=1;
Stoich(2,3)=-1;
%Rxn 3; Eh+A -> EAh, k2f
Stoich(3,3)=-1;
Stoich(3,4)=1;
%Rxn 4; EAh -> EA+h, k2r
Stoich(4,3)=1;
Stoich(4,4)=-1;
%Rxn 5; EAh -> Em + Q, k
Stoich(5,4)=-1;
Stoich(5,5)=1;

%Initialize the starting state with some hemimethylated sites,
%according to fracmeth
SitesArray0=zeros(NumCpGSites,NumSpec);
InitRandArray=rand(NumCpGSites,1);
SitesArray0(InitRandArray<=fracmeth,2)=1;
SitesArray0(:,1)=abs(SitesArray0(:,2)-1);
CurrState=SitesArray0;

[a,a_0]=GetPropensity(CurrState);

Time=0;
timeind=1;
Flag=0;
count=0;
%rng('shuffle')
while Time<TimePt
    %Gillespie algorithm: 2 random numbers, first to get time to next
    %reaction, 2nd to select which reaction
    Rands_KMC=rand(1,2);
    tau_KMC=1/a_0*log(1/Rands_KMC(1)); %time to next reaction
    [tau_FPT,sID,nID]=FPTStep(CurrState,sites);
    
    %the next event (KMC or FPT) corresponds to min tau
    if tau_KMC<=tau_FPT
        %disp('KMC')
        tau=tau_KMC;
        count=count+1;%count the number of steps
        if isinf(tau) %if no reaction is possible, increment the time anyway
            Time=Time+.5;
        else
            Time=Time+tau;
            findnextrxn=cumsum(a(:))>Rands_KMC(2)*a_0;
            NextRxn=sum(1-findnextrxn)+1; %index of next reaction
            [I,J]=ind2sub(size(a),NextRxn); %I returns the row index where the next
            %reaction occurs, J is the index of next reaction from NumRxns
            
            CurrState(I,:)=CurrState(I,:)+Stoich(J,:);
        end
    else
        %disp('FPT')
        tau=tau_FPT;
        count=count+1;
        if isinf(tau) %if no reaction is possible, increment the time anyway
            Time=Time+.5;
        else
            Time=Time+tau;
            %change the state of the start and stop sites for hop events
            CurrState(sID,5)= CurrState(sID,5)-1;
            CurrState(sID,6)= CurrState(sID,6)+1;
            if nID>0
                %if nID=0, then the enzyme goes back to solution, no change
                %to neighbor state
                CurrState(nID,2)= CurrState(nID,2)-1;
                CurrState(nID,3)= CurrState(nID,3)+1;
            end
        end
    end
    
    if TrajFlag==1
        savestep=[Time,sum(CurrState,1)];
        Traj(count,:)=savestep;
        %   AllTraj(:,:,count)=CurrState; %for checking processive mech, it may be
        %desired to save the entire per-site trajectory information
    end
    
    %Reset the propensity array and NeighborArray according to new state of
    %system
    [a,a_0]=GetPropensity(CurrState); %call the function that computes the propensity for each site
    %end
end
Time
SimDat(:,1)=SimDat(:,1)+CurrState(:,5)+CurrState(:,6);
SimDat(:,2)=SimDat(:,2)+CurrState(:,1)+CurrState(:,2)+CurrState(:,3)+CurrState(:,4);


    function [a,a_0]=GetPropensity(CurrState)
        %Calc the per-site propensities for each reaction. Array is
        %DNASize x NumRxns.
        Enum=E-sum(CurrState(:,3))-sum(CurrState(:,4))-sum(CurrState(:,5)); %number of free enzyme copies
        Anum=A-sum(CurrState(:,4));% why these others?-sum(CurrState(:,5))-sum(CurrState(:,6)); %number of AdoMet copies
        %Calc the per-site propensities for each reaction.
        a=zeros(NumCpGSites,NumRxns);
        a(:,1)=k1f*Enum*CurrState(:,2); %E + h -> [Eh]
        a(:,2)=k1r*CurrState(:,3); %  [Eh] -> E + h
        a(:,3)=k2f*CurrState(:,3)*Anum; %[Eh]+ A -> [hEA]
        a(:,4)=k2r*CurrState(:,4); %[hEA] -> [EA]+ h
        a(:,5)=kcat*CurrState(:,4); %[hEA] -> Em + Q
        
        a_0=sum(a(:));
    end

    function [tau_FPT,sID,nID]=FPTStep(CurrState,sites)
        %returns the sampled tau_FPT for super hops, and sID (the CpG ID of
        %the hop start, and nID (the CpG ID of the hop end)
        
        %first need a list of triples: all Em Site IDs and their nearest
        %neighbors in h-state to left and right (hopping targets)
        Curr_Em_inds=find(CurrState(:,5));
        %if there are no sites in Em state, the FPT step is moot, set FPT
        %to infinity so the event doesn't happen
        if isempty(Curr_Em_inds)
            tau_FPT=Inf;
            sID=0;
            nID=0;
        else
            Curr_Em_sites=sites(Curr_Em_inds);%sites(CurrState(:,5));
            Curr_h_inds=find(CurrState(:,2));
            Curr_h_sites=sites(Curr_h_inds);
            
            N_Em=numel(Curr_Em_sites);
            repC=repmat(Curr_h_sites(:)',N_Em,1);
            ga=repC<Curr_Em_sites(:);
            [Neighbor_L,track]=max(repC.*ga,[],2);
            Neighbor_Li=Curr_h_inds(track);
            ga=repC>Curr_Em_sites(:);
            za=repC.*ga;
            MaxSiteID=max(sites); %used later for vectorization
            za(za==0)=MaxSiteID+MaxDist+1; %trick used for vectorization
            [Neighbor_R,track]=min(za,[],2);
            Neighbor_Ri=Curr_h_inds(track);
            DiffSites_a=[Curr_Em_sites(:),Neighbor_L,Neighbor_R];
            %need to also keep track of indices in simulation array
            %DiffSites_i=[Curr_em_inds,Neighbor_Li,Neighbor_Ri];
            
            %make sure to only keep sites who have an available hop < MaxDist
            DistL=abs(Curr_Em_sites-Neighbor_L);
            DistR=abs(Curr_Em_sites-Neighbor_R);
            %%%
            KeepInds=find(sum([DistL<MaxDist,DistR<MaxDist],2));
            
            if isempty(KeepInds)
                %if there are no sites that are reachable by hopping, set the
                %FPT as exponential with rate koff (weak diffusion limit)
                tau_FPT=1/koff*log(1/rand); %sample from exponential (just like KMC)
                %randomly select an Em that it unbinds from
                sID=Curr_Em_inds(randi(length(Curr_Em_inds)));
                nID=0;
            else
                
                %set up domain info (L=size of domain, de=distance from left edge)
                Domains=[DiffSites_a(KeepInds,3)-DiffSites_a(KeepInds,2),DiffSites_a(KeepInds,1)-DiffSites_a(KeepInds,2)];
                N_Dom=numel(KeepInds);
                DiffSites=DiffSites_a(KeepInds,:);
                DiffSites_i=[Curr_Em_inds(KeepInds,:),Neighbor_Li(KeepInds,:),Neighbor_Ri(KeepInds,:)];
                
                %%%%
                %flip them if the current Em is closer to right edge (for looking up FPT
                %function)
                Domains_adj=[Domains(:,1),min(Domains(:,2),(Domains(:,1)-Domains(:,2)))];
                %Record whether L/R neighbor is closer: if CloseL==0 (right is closer), will need
                %to flip the ordering again later
                CloseL=Domains(:,2)<=Domains(:,1)-Domains(:,2);
                %Now Sample from the CDFs for to get exit time. Also sample
                %from Exit_Probs to get exit location (does it hop L, R, or
                %unbind back to solution with koff)
                Rands=rand(N_Dom,2);
                for ii=1:N_Dom
                    inds=Domains_adj(ii,:);
                    %this is an approximation: if the domain size is too large>MaxL=MaxDist*2,
                    %but we have kept it because the enzyme is <MaxDist from one edge, just assume the whole domain is MaxL
                    inds(1)=min(inds(1),MaxDist*2);
                    %CDF=CDFs{inds(1),inds(2)};
                    %this needs to be adjusted if not reading from file.
                    CDF=CDFs.CDFs{inds(1),inds(2)};
                    %ts=FPT_ts{inds(1),inds(2)};
                    ts=FPT_ts.FPT_ts{inds(1),inds(2)};
                    %sampling the time
                    get_ts=ts(CDF<=Rands(ii,1));
                    if isempty(get_ts)
                        sampt(ii)=eps;
                        t_ind=1;
                    else
                        sampt(ii)=get_ts(end); %array of FPTs
                        t_ind=numel(get_ts);
                    end
                    
                    %get the uncorrected Exit Probs
                    UnExitP=Exit_Probs.Exit_Probs{inds(1),inds(2)}(:,t_ind);
                    %correct the ExitPs: flip the Exit Probs if R neighbor is closer
                    ExitP=[UnExitP(2-CloseL(ii)),UnExitP(1+CloseL(ii)),UnExitP(3)];
                    %get corresponding exit states based on Exit Probs
                    findexitstate=cumsum(ExitP(:))>Rands(ii,2);
                    %exitstate =1 (hopped left) =2 (hopped R) =3 (unbind)
                    exitstate(ii)=sum(1-findexitstate)+1;
                end
                
                %get the minimum FPT
                [tau_FPT,samp_ind]=min(sampt);
                
                %which was the site of the Em that underwent diffusion?
                %ChosenSite=DiffSites(samp_ind,1);
                sID=DiffSites_i(samp_ind,1);
                %where did it hop to?
                %HopSite=DiffSites(samp_ind,3-ExitLeft(samp_ind));
                if exitstate(samp_ind)>2
                    %set to 0 if it went back to solution
                    nID=0;
                else
                    %get the neighbor index if it hopped L or R
                    nID=DiffSites_i(samp_ind,exitstate(samp_ind)+1);
                end
            end
        end
        
    end
end





