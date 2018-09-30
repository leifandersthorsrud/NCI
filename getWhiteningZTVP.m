function zstar=getWhiteningZTVP(z,beta,numDynFac,numVarTrans)

if size(beta,2)>1
    error([mfilename ':input'],'beta has more than one lag. Not supported')
end

zlagadj=repmat(beta,[1 size(z,2) size(z,3)]).*z;
zstar=z(:,1:numVarTrans,2:end);
cnt=numVarTrans;
for i=1:numDynFac                
    zstar=cat(2,zstar,z(:,1+cnt:numVarTrans+cnt,2:end)-zlagadj(:,1+cnt-numVarTrans:cnt,1:end-1));
    cnt=cnt+numVarTrans;                
end
zstar=cat(2,zstar,-zlagadj(:,1+cnt-numVarTrans:cnt,1:end-1));                                    