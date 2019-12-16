function [CDFs,FPT_ts,Exit_Probs] = MakeFPTLookup(D,koff,h,MaxDist)

MaxL=MaxDist*2;
%initialize the cell arrays to store all the computed CDFs
CDFs=cell(MaxL,MaxDist);
FPT_ts=cell(MaxL,MaxDist);
Exit_Probs=cell(MaxL,MaxDist);

%set up a tridiagonal matrix of max size so we don't have to keep building
%it each time for different domain sizes
Rmax=gallery('tridiag',MaxL+1,D/h^2,-2*D/h^2-koff,D/h^2);

for ii=2:MaxL
    MaxdeL=ceil(ii/2);
    L=ii
    Sz=L+2; %this makes the matrix the right size to have two targets (L and R) plus
    %an absorbing state which is the solution
    R=zeros(Sz);
    %pull the rate matrix out of Rmax with correct size
    R(1:Sz-1,1:Sz-1)=Rmax(1:Sz-1,1:Sz-1);
    
    %put in the unbinding rates
    R(Sz,:)=koff;
    %correct diagonal entries
    R=R-diag(diag(R));
    R=R-diag(sum(R,1));
    
    %make the ends and solution absorbing
    R(:,Sz)=0; %solution is absorbing
    R(:,1)=0;
    R(:,Sz-1)=0;
    
    inds_I=2:L;
    K_sub=R(inds_I,inds_I)';
    RHS=ones(numel(inds_I),1)*-1;
    % %TxB is the mean time to reach B from any x not in B
    TxB=linsolve(K_sub,RHS); % solving the linear system to get TxB
    
    for jj=1:MaxdeL
        MFPT_est=TxB(jj);
        t_end=10*MFPT_est; %number of times the estimated MFPT to go to
        t_mid=2*MFPT_est;
        ntimes=3000;
       % ntimes=5000;
        
        %subsume all absorbing states into first position for FPT calc
        R_one=R;
        R_one(1,:)=R_one(1,:)+R_one(Sz-1,:)+R_one(Sz,:);
        R_one=R_one(1:Sz-2,1:Sz-2);
        R_one(:,1)=0;
        %the absorbing state (b) is now always in position 1
        b=1;
        a=jj+1; %the matrix index of the starting site
t_early=linspace(0,t_mid,ntimes/2);       
t_axis=[t_early(1:end-1),linspace(t_early(end),t_end,ntimes/2)];


  [F,ts]=CalcFPT_RateMatrix(R_one,a,b,t_axis);

      %  [F,ts]=CalcFPT_RateMatrix(R_one,a,b,t_end,ntimes);
        CDFs{ii,jj}=cumsum(F);
        FPT_ts{ii,jj}=ts;
        
        %simulate
        P0=zeros(size(R,1),1);
        P0(a+1)=1;
        
        odefun=@(t,y) R*y;
        [t,y]=ode15s(odefun,t_axis,P0);
        %record the probability of exit L, R, or to solution
       %this is measured as relative probability flux to each exit state as fcn of time
       dP=R*y';
       Fluxes=dP([1,Sz-1,Sz],:);
       sumFluxes=repmat(sum(Fluxes,1),3,1);
       RelFluxes=Fluxes./sumFluxes; 
        Exit_Probs{ii,jj}=RelFluxes;%[y(end,1),y(end,Sz-1),y(end,Sz)];      
        
    end
end

end
