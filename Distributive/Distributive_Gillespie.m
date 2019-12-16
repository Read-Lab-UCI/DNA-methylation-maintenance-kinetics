function [SimDat] = Distributive_Gillespie(sites,fracmeth,TimePts)
%clear
tic
NumCpGSites=numel(fracmeth);
%the reaction scheme is E + A <-> [EA]+ h <-> [hEA] -> m + E + A', rates
%k1f,k1r,k2f,k2r,k3. The model is based on a modified compulsory-order
%ternary-complex mechanism in which A binds before h
%where h= hemimethylated substrate, m= methylated substrate
%E= DNMT1 enzyme, [EA]= enzyme-substrate complex 1
%A= AdoMet, [hEA]=enzyme-substrate-AdoMet complex 2
%k1f - forward rate, complex 1 formation
%k1r - reverse rate, complex 1 de-formation
%k2f - forward rate, complex 2 formation
%k2r - reverse rate, complex 2 de-formation
%k3 - catalytic step (product formation)

%KmA=k3/k1f
%Kmh=(k2r+k3)/k2f
%kcat=k3
%k1r=k1f/10;%arbitrary value

%Enzymatic parameters (KmA, Kmh and kcat) of the model are taken from Pradhan et al. (1999). DNMT1 with SNRPN exon-1 as substrate.
%K are given in micromolar



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

% k1f=kcat/Kmh; % /(uM*h)
% k2f=(kcat)/(KmA-1/10); %/(uM*h)
% k2r=k2f/10; %arbitrary relationship  (/hr)
% k1r=k1f/10; %arbitrary relationship   (/hr)


k1f=kcat/Kmh;
k2r=100; %arbitrary value(/hr)
k1r=Kia*k1f;
k2f=(kcat+k2r)/(KmA);
k3=kcat;

%following calculations are based on assumption of nuclear volume=1000
%microns-cubed=1E-12 Liters:
V_nuc=1E-12;
E=round(E0*1E-6*6.022E23*NumCpGSites*V_nuc/h); % # of enzyme copies scaled by NumCpGSites
A=round(A0*1E-6*6.022E23*NumCpGSites*V_nuc/h); % # of AdoMet copies scaled by NumCpGSites

% Kmh=Kmh*1E-6*6.022E23*V_nuc;% # of copies scaled by NumCpGSites
% KmA=KmA*1E-6*6.022E23*V_nuc;% # of copies scaled by NumCpGSites
k1f=k1f*h/NumCpGSites/V_nuc*1E6/6.02E23; % /(copies*h) scaled by NumCpGSites
k2f=k2f*h/NumCpGSites/V_nuc*1E6/6.02E23; % /(copies*h) scaled by NumCpGSites


NumRxns=5;
NumSpec=4; % [u,h,EAh,m] EA is considered differently since it is not a possible state of a CpG

SimDat=zeros(NumCpGSites,2); %initialize AllDat, 2nd dimension is # methylated, #
%unmethylated

Stoich_CpG=zeros(NumRxns,NumSpec); %stoichiometry of reactions in terms of [u,h,EAh,m]

Stoich_CpG(3,2)=-1;
Stoich_CpG(3,3)=1;
Stoich_CpG(4,2)=1;
Stoich_CpG(4,3)=-1;
Stoich_CpG(5,3)=-1;
Stoich_CpG(5,4)=1;

Stoich_EA=zeros(NumRxns,1); %stoichiometry of reactions in terms of EA
Stoich_EA(1)=1;
Stoich_EA(2)=-1;
Stoich_EA(3)=-1;
Stoich_EA(4)=1;



%Initialize the starting state with some hemimethylated sites
SitesArray0=zeros(NumCpGSites,NumSpec); %states [u,h,EAh,m,A]
InitRandArray=rand(NumCpGSites,1);
AssignedFrac=InitRandArray<=fracmeth;
SitesArray0(AssignedFrac,2)=1; %identify the sites that are h in the starting state
SitesArray0(:,1)=abs(SitesArray0(:,2)-1); %the sites that are not h in the starting state are
EA=0;


CurrState=SitesArray0;
Enum=E-EA-sum(CurrState(:,3)); %number of free enzyme copies
Anum=A-EA-sum(CurrState(:,3)); %number of AdoMet copies
%Calc the per-site propensities for reactions 3 to 5
a_3to5=zeros(NumCpGSites,NumRxns-2);
a_3to5(:,1)=k2f*CurrState(:,2)*EA; %[EA]+ h -> [hEA]
a_3to5(:,2)=k2r*CurrState(:,3); %[hEA] -> [EA]+ h
a_3to5(:,3)=k3*CurrState(:,3); %[hEA] -> m + E + A'

%Calc the per-site propensities for reactions 3 to 5
a_1to2=zeros(2,1);
a_1to2(1)=k1f*Enum*Anum; %E + A -> [EA]
a_1to2(2)=k1r*EA; %  [EA] -> E + A

a_0=sum(a_1to2)+sum(a_3to5(:)); %Total sum of promensities

Time=0;
timeind=1;
count=0;
Flag=0;

while Time<TimePts(end)
    Rands=rand(1,2); %pick 2 random numbers
    tau=1/a_0*log(1/Rands(1)); %the first number is used to set the time-step following Gillespie algorithm
    count=count+1; %count the number of steps
    %a_0_count(count,1)=a_0;
    if isinf(tau)
        Time=Time+.5;
    else
        Time=Time+tau;
        findnextrxn=cumsum([a_3to5(:);a_1to2])>Rands(2)*a_0; %reactions 1 and 2 are placed at the bottom of the cumsum array
        NextRxn=sum(1-findnextrxn)+1; %find the position in the cumsum array of the reaction that takes place
        
        [I,J]=ind2sub(size(a_3to5),NextRxn); %find the position in a of the reaction that takes place (J) and the site in which takes places (I)
        
        
        if J<4 % If the reaction that takes place is 3, 4 or 5, J=<3 (Nrows of a_3to5)
            J=J+2;
            CurrState(I,:)=CurrState(I,:)+Stoich_CpG(J,:); %reactions 3,4 and 5 modify both the state of the CpGs and EA
            EA=EA+Stoich_EA(J);
            
        else % If the reaction that takes place is 1 or 2, J=4 (Nrows of a_3to5)
            if I==1 % If the reaction that takes place is 1, I=1
                J=J-3;
            else  % If the reaction that takes place is 2, I=2
                J=J-2;
            end
            EA=EA+Stoich_EA(J); %reactions 1 and 2 only modify the state of EA
        end
        
    end  
    
    SimDat(:,1)=SimDat(:,1)+CurrState(:,4); %methylated reads
    SimDat(:,2)=SimDat(:,2)+CurrState(:,1)+CurrState(:,2)+CurrState(:,3); %unmethylated reads
    
    
    Enum=E-EA-sum(CurrState(:,3)); %Recalculate the number of free enzyme copies
    Anum=A-EA-sum(CurrState(:,3)); %Recalculate the number of AdoMet copies
    
    
    %     Recalc the per-site propensities for each reaction
    a_3to5(:,1)=k2f*CurrState(:,2)*EA; %[EA]+ h -> [hEA]
    a_3to5(:,2)=k2r*CurrState(:,3); %[hEA] -> [EA]+ h
    a_3to5(:,3)=k3*CurrState(:,3); %[hEA] -> m + E + A'
    
    a_1to2=zeros(2,1);
    a_1to2(1)=k1f*Enum*Anum; %E + A -> [EA]
    a_1to2(2)=k1r*EA; %  [EA] -> E + A
    
    a_0=sum(a_1to2)+sum(a_3to5(:));
end

end


%save(datafilename,'AllDat','sites');
%toc
%display (Traj(:,:,1))
%display (a_0_count)