function z=getAasCoefficient(nvarObs,A,stacked)

[nobs,numberOfStates]=size(A);

if stacked
    z=repmat(permute(A',[3 1 2]),[nvarObs 1 1]);
else
    z=zeros(nvarObs,numberOfStates*nvarObs,nobs);
    cnt=1;
    for i=1:nvarObs
        z(i,cnt:cnt+numberOfStates-1,:)=A';
        cnt=cnt+numberOfStates;
    end     
end