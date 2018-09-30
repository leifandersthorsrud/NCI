function [P, H, Z, ZU, W, D] = bmffvar_doObsEst(y, AirepLargeL, P, Z, D, H, W,...
                priorSettings, numDynFac, freqIndex, restrict, freqType)                    
                    
scale       = 1;                        
w           = priorSettings.W(freqIndex, freqIndex);
vw          = ceil(priorSettings.vW*scale);
d           = priorSettings.d(freqIndex);
vh          = ceil(priorSettings.vH*scale);
h           = priorSettings.H(freqIndex, freqIndex);
p           = priorSettings.P(freqIndex, freqIndex);
vp          = priorSettings.VP(freqIndex, freqIndex);

[P, Z, ZU, D, H, W] = ...
               doObsEstFreqSpec(y, P, Z, H, W, D, AirepLargeL,...
               numDynFac, w, vw, d, vh, h, p, vp, restrict, freqType);                  
                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, Z, ZU, D, H, W] = ...
               doObsEstFreqSpec(y, P, Z, H, W, D, AL,...
               numDynFac, w, vw, d, vh, h, p, vp, restrict, freqType)                  
           
[nObs, N]   = size(y);
nlagAr      = size(P, 2);
numDynFacZ  = numDynFac + 1;

yw  = getWhiteningY( y, P);
x   = permute(getWhiteningA( AL, P, numDynFac, 1 ), [3 2 1]);

H   = diag(H);
W   = reshape(diag(W), [numDynFacZ N]);
ZU  = Z;

w   = reshape(diag(w), [numDynFacZ N]);
h   = diag(h);
d   = reshape(d, [numDynFacZ N]);
p   = diag(p);
vp  = diag(vp);

parfor i = 1 : N        
    
    % 1) Pass the variables
    bdraw               = permute(Z(i, :, :), [3 2 1]);    
    ddraw               = D(i);        
    hdraw               = H(i);    
    
    vqi                 = vw;       
    qi                  = w(i);   
    di                  = d(i);                      
    hi                  = h(i);
    pi                  = p(i);
    vpi                 = vp(i);
    
    ywi                 = yw(:, i);    
    xi                  = x(:, :, i);
        
    nonanIdx            = ~isnan(ywi);        
    
    if i == 1 && restrict 
        
        wdraw = 0;
        ddraw = 0;
        % Note: using nan here. This will be filled in below
        [bdraw, bdrawu] = deal(nan(nObs, 1)); %ones
        
        if strcmpi(freqType, 'd')
            % 1) Restrict daily factor loading to 1 for first element            
                        
            e = ywi(nonanIdx) - xi(nonanIdx, 1 : numDynFacZ);
            
            % Compute error
            v           = nObs + vh;            
            tohdraw     = (e'*e + vh*hi)/v;
            hdraw       = 1/gamm_rnd(1, 1, .5*v, .5*v*tohdraw);                            
            
            bdraw(:, 1)     = 1;
            bdrawu(:, 1)    = 1;
                                    
        else
            % 2) Restrict first element among non-daily variables to be a
            % constant (no-tvp)                    
            [hdraw, bdrawi]     = bmffvar_drawB(...
                                    ywi(nonanIdx), xi(nonanIdx, 1 : numDynFacZ),...
                                    hdraw, vh, hi,...
                                    1, 0.1);            
            bdraw(:, 1)  = bdrawi;
            bdrawu(:, 1) = bdrawi;                                    
            
        end
        
    else                
        
        if strcmpi(freqType, 'd')
            % 3) Do tvp for daily factor loadings 
            [bdrawi, bdrawui, wdraw, hdraw, ddraw,...
                ~, ~] = dynamicLatentThresholdBayes(ywi(nonanIdx), xi(nonanIdx, 1 : numDynFacZ), bdraw(nonanIdx, :, :),...
                                 vqi, qi,...
                                 ddraw, di,...
                                 hdraw, vh, hi);                                                        
            
        else                                    
            % 4) Do constant loading for other frequencies. Need this due
            % to the aggregation restrictions. Done here in separate step
            % relative to 2) above because of the different prior
            % assumption. 
            [hdraw, bdrawi] = bmffvar_drawB(...
                                    ywi(nonanIdx), xi(nonanIdx, 1 : numDynFacZ),...
                                    hdraw, vh, hi,...
                                    0, 1);
            
            bdrawi  = repmat(bdrawi, [sum(nonanIdx) 1]);
            bdrawui = bdrawi;
            
            wdraw   = 0;
            ddraw   = 0;           
        end
                             
        bdraw(nonanIdx, :)  = bdrawi;
        bdrawu              = zeros(nObs, numDynFacZ);
        bdrawu(nonanIdx, :) = bdrawui;
                
    end        
    
    Z(i, :, :)   = bdraw';
    ZU(i, :, :)  = bdrawu';
    W(:, i)      = wdraw;
    H(i)         = hdraw;
    D(i)         = ddraw;
             
    
    if any(strcmpi(freqType, {'q', 'm', 'w'})) && i == 1            
        P(i, :) = 0;        
    else        
        nonanIdx            = ~isnan(y(:, i));      
        ALi                 = AL(nlagAr + 1 : end, 1);
        Pi      = doARdraw(pi, vpi, y(nonanIdx, i), ALi(nonanIdx), bdrawu(nonanIdx, :), nlagAr, hdraw);        
        P(i, :) = Pi;                            
        
    end
    
end

% Adjust output again
H = diag(H);
W = diag(vec(W));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = doARdraw(p, vp, y, x, b, nlag, hdraw)
        

capvp0inv   = inv(vp);
eraw        = y - x.*b;           
e           = eraw(nlag + 1 : end, 1);

bige        = latMlag(eraw, nlag);
bige        = bige(nlag + 1 : end, :);
bigesquare  = bige'*bige;

capvp1inv   = capvp0inv + hdraw*bigesquare;
capvp1      = inv(capvp1inv);
bar1        = capvp1*(capvp0inv*p' + hdraw*bige'*e); 
%Check and see if draw satisfies stationarity, if not redraw
stabilityCondition = 0;
while stabilityCondition == 0
    bardraw             = bar1 + norm_rnd(capvp1);            
    [bardrawComp, ~]    = varGetCompForm(bardraw', [], nlag, 1);                        
    if any(abs(eig(bardrawComp)) >=1 )
        stabilityCondition = 0;
    else
        stabilityCondition = 1;
    end
end

P = bardraw;