clear
NumCpGSites=100000;
load('WTArrest_Bulk_Chr1_Processed_curated.mat') %load methylation data for initialization
start=1;
CpGFileInds=start:start+NumCpGSites-1;
Chunk=10000;
AllDat=zeros(Chunk,4,2); %initialize simulated data array, 4 exp timepoints
NChunks=NumCpGSites/Chunk;

Sites=Sites(CpGFileInds);
Frac=Frac(CpGFileInds);
for loopchunk=1:NChunks
    Inds=Chunk*(loopchunk-1)+1:Chunk*loopchunk;
    sites=Sites(Inds);
    fracmeth=Frac(Inds);
    save('NumCpGSites', 'NumCpGSites');
    
    Times=[0.5,1.5,4.5,16.5]; %Requested time-points in hours
    NumReads_range=[10,5,5,5]; %number of reads requested at each of associated timepoints
    TrajFlag=0;
    SimFlag=1;
AllDat=zeros(Chunk,4,2); %initialize simulated data array, 4 exp timepoints
    for j=1:numel(Times)
        TimePts=Times(j);
        NumReads=NumReads_range(j);
        datafilename=['AllDat_'  'Simulated_FPTKMC' num2str(NumCpGSites) '_Ch' num2str(loopchunk) ];
        
        for loopreads=1:NumReads
            [AllDat_each]=Processive_FPTKMC_1DDiff(sites,fracmeth,TimePts,TrajFlag,SimFlag);
            
            if loopreads==1
                AllDat_time=AllDat_each;
            else
                AllDat_time=AllDat_time+AllDat_each;
            end
            SimFlag=1;
        end
        
        AllDat(:,j,:)=AllDat(:,j,:)+reshape(AllDat_time,[Chunk,1,2]);
        
    end
    save(datafilename,'AllDat','sites');
end
