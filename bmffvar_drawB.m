function [Qirep, bdraw] = bmffvar_drawB(e, a, Qirep, vprior, Qprior, T, V)

if nargin == 5
    V = 0.1;
    T = 1;  
end


betan       = a\e;
hdraw       = 1/Qirep;
capv0inv    = inv(V(1, 1));
xsquare     = a'*a;
v1          = vprior + size(e, 1);    
v0s02       = vprior*Qprior;
b0          = T;
%draw from beta conditional on h
capv1inv    = capv0inv + hdraw*xsquare;
capv1       = inv(capv1inv);
b1          = capv1*(capv0inv*b0 + hdraw*xsquare*betan);
bdraw       = b1 + norm_rnd(capv1);       
%draw from h conditional on beta
s12         = ((e - a*bdraw)'*(e - a*bdraw) + v0s02)/v1;
hdraw       = gamm_rnd(1, 1, .5*v1, .5*v1*s12);    
Qirep       = 1/hdraw;

