function [Q, T, BiQ] = bmffvar_doTransEst(nlag, TRes, Airep, Qirep,... 
                            BiQ, priorSettings)
                                                
nobs        = size(Airep, 1) - nlag;
constant    = false;

[Y, X, K]   = getYZ(Airep, nlag, constant, TRes);        

n       = 1;
H       = sparse(nobs*n, nobs*n);
idx     = 1:nobs:nobs*n - n;
for i = 1 : nobs        
    ht = inv(Qirep(:, :, i));                        
    for k = 1 : n            
        H(idx, idx) = ht;        
    end
    idx = idx + 1;
end;

[T, ~]    = transitionSimBetaBayes(nobs, nlag, constant, ...
                        K, Y, X, H, priorSettings.VT', priorSettings.T', TRes);                                        

%% Draw error term

ehat    = Airep - latMlag(Airep, nlag)*T(1, :)';                    
ehat    = ehat(nlag + 1 : end, :);    
Sigmai  = permute(sqrt(Qirep), [3 1 2]);

a0Q = log(priorSettings.a0Q);
[Sigmai, BiQ, ~] = drawSigmaTVP(ehat, Sigmai, BiQ,...
                        a0Q, priorSettings.p0Q, priorSettings.vB, priorSettings.B);                   
Q               = Sigmai.^2;
    
