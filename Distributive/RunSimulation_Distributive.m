clear
NumCpGSites=10000;
load('WTArrest_Bulk_Chr1_Processed_curated.mat') %load methylation data for initialization
start=1;
CpGFileInds=start:start+NumCpGSites-1;
Chunk=5000;
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
    AllDat=zeros(Chunk,4,2); %initialize simulated data array, 4 exp timepoints
    for j=1:numel(Times)
        TimePts=Times(j);
        NumReads=NumReads_range(j);
        datafilename=['AllDat_'  'Simulated_Distributive' num2str(NumCpGSites) '_Ch' num2str(loopchunk) ];
        
        for loopreads=1:NumReads
            [AllDat_each]=Distributive_Gillespie(sites,fracmeth,TimePts);
            
            if loopreads==1
                AllDat_time=AllDat_each;
            else
                AllDat_time=AllDat_time+AllDat_each;
            end
        end
        
        AllDat(:,j,:)=AllDat(:,j,:)+reshape(AllDat_time,[Chunk,1,2]);
        
    end
    save(datafilename,'AllDat','sites');
end
