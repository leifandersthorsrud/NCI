function zstar=getWhiteningA(A,beta,numDynFac,numVarTrans)

[nvarObs,nlag]=size(beta);
[nobs,numberOfStates]=size(A);

if nlag>1
    error([mfilename ':input'],'beta has more than one lag. Not supported')
end;

% Get A as coefficient matrix (nvarObs x numberOfStates x t)
z=getAasCoefficient(nvarObs,A,true);

zlagadj=repmat(beta,[1 numberOfStates nobs]).*z;
zstar=z(:,:,2:end)-zlagadj(:,:,1:end-1);


