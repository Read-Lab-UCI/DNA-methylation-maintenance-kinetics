function [S,ts] = CalcFPT_RateMatrix(R,a,b,ts)%t_end,ntimes)
%input arguments:
%R the rate matrix in proper format (columns sum to 0)
%a the index of the initial state
%b the index of the final state
%t_end the end time of the distribution
%ntimes the number of computed timepoints

%ts=linspace(0,t_end,ntimes);
M=-R; %take negative to be consistent with Gillespie
M(:,b)=0; %make the b state absorbing--all column entries =0

[V,DD]=eig(M);
lambda=diag(DD);
Vi=inv(V);

getinds=find(lambda>0);

C=-V(b,:).*Vi(:,a)';
ntimes=length(ts);
S=zeros(length(getinds),ntimes);
expmat=exp(-lambda(getinds)*ts);
lamc=lambda(getinds).*C(getinds)';
Rlamc=repmat(lamc,1,ntimes);
Smat=Rlamc.*expmat;
S=sum(Smat,1);
S=S/sum(S);

end
