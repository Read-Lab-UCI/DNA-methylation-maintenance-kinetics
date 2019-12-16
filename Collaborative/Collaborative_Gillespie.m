function [SimDat] = Collaborative_Gillespie(sites,fracmeth,TimePt,TrajFlag)

%input arguments:
%AssignedFrac is an Nx1 array with values in [0,1], where
%N=Number of CpG sites to be simulated in the model and values are the
%assigned probability of steady state methylation at that site.

%Sites is an Nx1 array with the basepair index of each site to be simulated

%TimePt is the requested TimePt for state of system (in hours)

%TrajFlag. If 1, the full trajectory will be saved. If 0, only the AllDat
%file with state for each TimePts will be returned.

NumCpGSites=numel(fracmeth); %number of CpG sites to be simulated
nDist=1E9;
nN=200;
i_rec=16E3;

% NeighborDist=16; %need to define a genomic distance for the processive mechanism, i.e., if there is
% %a DNMT1 bound within this distance of another hemi-methylated site, it can
% %reach it processively. This is currently set arbitrarily (but short). This
% %is a discrete threshold mechanism, i.e., if the current site is within NeighborDist
% %of an enzyme-bound site, the enzyme can proceed to current site
% %processively, else not.
% save('NeighborDist','NeighborDist')

%COTPM adapted from Pradhan
%E + h <-> Eh k1f, k1r
%EA + h <-> EAh (I added these reactions) k2f, k2r

%Processive mechanism from Svedruzic + diffusion ideas from Erban
%EQm + h -> m + Eh     %DNMT1 hop to the right (k_hop_R)
%h + EQm -> Eh + m     %DNMT1 hop to the right (k_hop_L)

%EQmi-> koff off rates

%following calculations are based on assumption of nuclear volume=1000
%microns-cubed=1E-12 Liters, and [DNMT1]=5000 nM-500 microMolar.
V_nuc=1E-12;
h=28E6; % # CpG sites
Ratio=0.01;
h0=h/V_nuc/6.022E23*1E6;%initial [h] in micromolar
E0=h0*Ratio; %total number of DNMT1 enzymes in system
kcat=40; %/hr
Kmh=1.3;%uM
KmA=5.5;%uM
A0=100; %uM
Kih=2.85; %uM

% k1f=kcat/Kmh; % /(uM*h)
% k2f=(kcat)/(KmA-1/10); %/(uM*h)
% k2r=k2f/10; %arbitrary relationship  (/hr)
% k1r=k1f/10; %arbitrary relationship   (/hr)


k1f=kcat/Kmh;
k2r=100; %arbitrary value(/hr)
k1r=Kih*k1f;
k2f=(kcat+k2r)/(KmA);



E=round(E0*1E-6*6.022E23*NumCpGSites*V_nuc/h); % # of enzyme copies scaled by NumCpGSites
A=round(A0*1E-6*6.022E23*NumCpGSites*V_nuc/h); % # of AdoMet copies scaled by NumCpGSites
% Kmh=Kmh*1E-6*6.022E23*V_nuc;% # of copies scaled by NumCpGSites
% KmA=KmA*1E-6*6.022E23*V_nuc;% # of copies scaled by NumCpGSites
k1f=k1f*h/NumCpGSites/V_nuc*1E6/6.02E23; % /(copies*h) scaled by NumCpGSites
k2f=k2f*h/NumCpGSites/V_nuc*1E6/6.02E23; % /(copies*h) scaled by NumCpGSites

a=650/1E3; %Parameter extracted from Lovkist
b=196/1E3; %Parameter extracted from Lovkist

%D=10^9; %bp^2/h. Arbitrary value
%nN=100; %Max # of neighbors a site can have if nDist is below 3600 bp
%nN=1000; %Maximum distance the E can travel in a processive manner
%determine all the left neighboring distances of each site
k_rec_L=zeros(NumCpGSites,nN); %initialize an array to store the distance between
%each site and its right neighbor.
%The last CpG on the right has no neighbors at the right
k_rec_R=zeros(NumCpGSites,nN); %initialize an array to store the distance between
%each site and its left neighbor.
%The first CpG on the left has no left neighbor

%initializing distance arrays for collaborative hops
Mesh=meshgrid(sites);
DistanceMatrix=abs(Mesh-Mesh');
NeighborInds_R=zeros(NumCpGSites,nN);
NeighborInds_L=zeros(NumCpGSites,nN);
for ii=1:NumCpGSites
    RHS=min([NumCpGSites,ii+nN]);
    LHS=max([1,ii-nN]);
    elsL=ii+1:RHS;
    elsR=LHS:ii-1;
    k_rec_L(ii,1:numel(elsL))=(DistanceMatrix(ii,ii+1:RHS)+b).^-1*a*k1f*i_rec;
    k_rec_R(ii,end-numel(elsR)+1:end)=(DistanceMatrix(ii,LHS:ii-1)+b).^-1*a*k1f*i_rec;
    NeighborInds_L(ii,1:numel(elsL))=ii+1:RHS;
    NeighborInds_R(ii,end-numel(elsR)+1:end)=LHS:ii-1;
end
k_rec_R=fliplr(k_rec_R);
NeighborInds_L(NeighborInds_L==0)=NumCpGSites+1; %trick for vectorization
NeighborInds_R(NeighborInds_R==0)=NumCpGSites+1; %trick for vectorization

clear Mesh
clear DistanceMatrix
disp('finished distance calcs')
NumRxns=7;  %number of reactions
NumSpec=6; % [u,h,Eh,EAh,EQm,m] %species in the model

SimDat=zeros(NumCpGSites,2); %initialize AllDat, 2nd dimension is # methylated, #
%unmethylated

Stoich=zeros(NumRxns,NumSpec); %stoichiometry of reactions in terms of [u,h,Eh,EAh,EQm,m]

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

%Rxn 5; EAh -> E+Q +m, k
Stoich(5,4)=-1;
Stoich(5,6)=1;

%Rxn 6; h + EQm-> Eh + m  %DNMT1 hop to the left (k_hop_L).
Stoich(6,2)=-1;
Stoich(6,3)=1;


%Rxn 7; EQm + h -> m + Eh  %DNMT1 hop to the right (k_hop_R).
Stoich(7,2)=-1;
Stoich(7,3)=1;

%Initialize the starting state with some hemimethylated sites,
%according to fracmeth
SitesArray0=zeros(NumCpGSites,NumSpec); %states [u,h,EAh,m,A]
InitRandArray=rand(NumCpGSites,1);
SitesArray0(InitRandArray<=fracmeth,2)=1;
SitesArray0(:,1)=abs(SitesArray0(:,2)-1);
CurrState=SitesArray0;

[a,a_0]=GetPropensity(CurrState); %call the function that computes the propensity for each site

Time=0;
timeind=1;
Flag=0;
count=0;

while Time<TimePt
    %Gillespie algorithm: 2 random numbers, first to get time to next
    %reaction, 2nd to select which reaction
    Rands=rand(1,2);
    tau=1/a_0*log(1/Rands(1)); %time to next reaction
    count=count+1;%count the number of steps
    if isinf(tau) %if no reaction is possible, increment the time anyway
        Time=Time+.5;
    else
        Time=Time+tau;
        findnextrxn=cumsum(a(:))>Rands(2)*a_0; %reactions 1 and 2 are placed at the bottom of the cumsum array
        NextRxn=sum(1-findnextrxn)+1; %index of next reaction
        [I,J]=ind2sub(size(a),NextRxn); %I returns the site index where the next
        %reaction occurs, J is the index of next reaction from NumRxns
        
        if J<6
            CurrState(I,:)=CurrState(I,:)+Stoich(J,:);     
        end
        
        if and(6<=J,J<=(5+nN)) %the processive reaction: must take care of the leaving site
            %  If any og these reactios take place it means that the enzyme has
            %  jumped from I+(1:nN) to I (from right to left)
            CurrState(I,:)=CurrState(I,:)+Stoich(6,:);
            
            
            %                 Neighbor=I+J-5;
            %                 CurrState(Neighbor,:)=CurrState(Neighbor,:)+Stoich(8,:); %increment the stoichiometry
        end
        if and((6+nN)<=J,J<=(5+2*nN))
            %  If any of these reactions take place it means that the enzyme has
            %  jumped from I-(1:nN) to I (from left to right)
            CurrState(I,:)=CurrState(I,:)+Stoich(7,:);
            %                 Neighbor=I-(J-(5+nN));
            %                 CurrState(Neighbor,:)=CurrState(Neighbor,:)+Stoich(8,:); %increment the stoichiometry
            
        end
        
        
        
    end
    
    %         if Flag< 1 && Time>=TimePts(timeind)
    %             Time %display current TimePt
    %
    %             %populate the SimDat array with the current state of the system, in
    %             %terms of unmethylated on nascent strand (u,h,Eh states) and
    %             %methylated (Em and m).
    %             SimDat(:,timeind,1)=SimDat(:,timeind,1)+CurrState(:,5)+CurrState(:,6);
    %             SimDat(:,timeind,2)=SimDat(:,timeind,2)+CurrState(:,1)+CurrState(:,2)+CurrState(:,3)+CurrState(:,4);
    %             timeind=timeind+1;
    %             if Time>TimePts(end) %indicator to exit simulation if TimePts(end) exceeded
    %                 Flag=1;
    %             end
    %         end
    
    if TrajFlag==1
        savestep=[Time,sum(CurrState,1)];
        Traj(count,:)=savestep;
        %   AllTraj(:,:,count)=CurrState; %for checking processive mech, it may be
        %desired to save the entire per-site trajectory information
    end
    
    
    %Reset the propensity array and NeighborArray according to new state of
    %system
    [a,a_0]=GetPropensity(CurrState); %call the function that computes the propensity for each site
end
Time
SimDat(:,1)=SimDat(:,1)+CurrState(:,5)+CurrState(:,6);
SimDat(:,2)=SimDat(:,2)+CurrState(:,1)+CurrState(:,2)+CurrState(:,3)+CurrState(:,4);
%

% %save the trajectory file, if flag
% if TrajFlag==1
%     save Traj Traj
%   %  save AllTraj AllTraj
% end


    function [a,a_0]=GetPropensity(CurrState)
        %Calc the per-site propensities for each reaction. Array is
        %NumCpGSites x NumRxns.
        Enum=E-sum(CurrState(:,3))-sum(CurrState(:,4))-sum(CurrState(:,5)); %number of free enzyme copies
        Anum=A-sum(CurrState(:,4))-sum(CurrState(:,5))-sum(CurrState(:,6)); %number of AdoMet copies
        %Calc the per-site propensities for each reaction.
        %Calc the per-site propensities for reactions
        a=zeros(NumCpGSites,NumRxns+2*(nN-1));
        a(:,1)=k1f*Enum*CurrState(:,2); %E + h -> [Eh]
        a(:,2)=k1r*CurrState(:,3); %  [Eh] -> E + h
        a(:,3)=k2f*CurrState(:,3)*Anum; %[Eh]+ A -> [hEA]
        a(:,4)=k2r*CurrState(:,4); %[hEA] -> [EA]+ h
        a(:,5)=kcat*CurrState(:,4); %[hEA] -> E+Q+m
        
        %for processive mechanism, propensity depends on neighboring states
        %Hop from right to left
        NeighborState_L=zeros(NumCpGSites,nN);
        Dummy=CurrState(:,3);
        Dummy(end+1)=0;
        NeighborState_L=Dummy(NeighborInds_L);
        
        %NeighborState_L is an ixj array that in row i will display 1 if neighbor j on the right is Eh
        %i.e. it is a candidate for enzyme to jump from j to site i
        Prop_hop_L=NeighborState_L.*k_rec_L; %propensity of a site i to revieve E from any of it nN neighbors on the right
        a(:,6:5+nN)=CurrState(:,2).*Prop_hop_L*Enum;
        
        %Hop from left to right
        %NeighborState_R is an ixj array that in row i will display 1 if neighbor j on the left is Eh
        %i.e. it is a candidate for enzyme to jump from j to site i
        NeighborState_R=zeros(NumCpGSites,nN);
        NeighborState_R=fliplr(Dummy(NeighborInds_R));
        a(:,6+nN:5+2*nN)=CurrState(:,2).*NeighborState_R.*k_rec_R*Enum;
        
        a_0=sum(a(:));
    end
end



