function [bdraw, bdrawu, wdraw, hdraw, dd,...
    acceptb, acceptd] = dynamicLatentThresholdBayes(y, x, bdraw,...
                             vprior_w, Sprior_w, dd, ...
                             thresholdK, hdraw, vprior_h, Sprior_h)
% Replication codes for
% Nakajima, J. & West, M
% Bayesian Analysis of Latent Threshold Dynamic Models
% Journal of Business & Economic Statistics
% 2013, 31, 151-164
%
                         
                         
[nT, nN]    = size(x);
invh        = repmat(inv(hdraw),[1 1 nT]);

phi     = eye(nN);
mu      = zeros(nN, 1);

%% Do the drawings

% Conditional on bju, draw vdrawj
wdraw        = drawV(bdraw, vprior_w, Sprior_w);        

% Conditional on sigma and B_T, impose the MH step for B_T                        
[bdraw, acceptb]  = drawBUpdated(y, x, bdraw, invh, wdraw, phi, mu, dd);     

% Conditional above, draw threshold
vdc                   = diag(abs(mu)) +  diag(thresholdK); %thresholdK.*Sprior_w{i};% vprior_w{i}*Sprior_w{i}; 
[dd, acceptd]         = drawBThreshold(y, x, dd, bdraw, invh, vdc);            

% Construct bdraw where threshold is imposed
idx             = abs(bdraw) < repmat(dd',[nT 1]);             
vb              = vec(bdraw);
vb(vec(idx))    = 0;
bdrawu          = reshape(vb,[nT nN]);                                

% Conditional on the above, draw hdraw, i.e., the observaton equation
% variance
hdraw = drawH(y, x, bdrawu, vprior_h, Sprior_h);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdraw = drawH(y, x, b, vprior_h, Sprior_h)
    
nT      = size(y, 1);
emat    = zeros(nT, 1);   

for t = 1 : nT

    emat(t, 1) = y(t) - x(t, :)*b(t, :)';
    
end

v1          = vprior_h + nT;        
H           = inv(vprior_h*Sprior_h + emat'*emat);            
hdraw       = inv( wish_rnd(H, v1) );  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vdraw = drawV(B,vprior_v,Sprior_v)
    

nobs        = size(B,1);

emat        = B(2:end,:) - B(1:end - 1,:);
v1          = vprior_v + nobs;        
H           = inv(vprior_v*Sprior_v + emat'*emat);            
vdraw       = inv( wish_rnd(H, v1) );  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% /*
% **  fSampLTBv samples beta with latent threshold
% **  in multi-variate model by single-move sampler
% **
% **  [NOTE]
% **    y_t = X_t*b_t + e_t,    e_t ~ N(0, G2)
% **    b_t = B_t*s_t,          s_t = I(B_t>=d)
% **
% **    B_{t+1} = mu + Phi*(B_t - mu) + eta_t,
% **    eta_t ~ N(0, Sig),  eta_0 ~ N(0, Sig/(I-Phi^2))
% **
% */
% LTVAR::fSampLTBv(const my, const amX, const mPhi,
% 				 const vmu, const mSig, const amG2inv,
% 				 const vd, const mbo)
% {
% 	decl ns, nk, np, mba, mSigi, mSigp, vpm1, vpm2;
% 	decl mSig0i, mSighi, mxh, mSigh, vbh, vbn, vbo, vba_1;
% 	decl dhn, dho, dln, dlo, dfrac, ca, i;

    
function [mba, ca] = drawBUpdated(my,amX,mbo,amG2inv,mSig,mPhi,vmu,vd)    
    
%ns = rows(my);			//# of sample
%nk = columns(my);		//# of series
%np = columns(amX[0]);	//# of states

ns  = size(my,1);    	
nk  = size(amX,2);
np  = nk;	
% number of accepted draws
ca  = 0;

% set porposal state equal to current state vector
mba = mbo;

mSigi   = inv(mSig);		% //Sig^(-1)
mSigp   = mSigi * (eye(np) + mPhi.^2);
vpm1    = (eye(np) - mPhi) * vmu;
vpm2    = (eye(np) - 2 * mPhi + mPhi.^2) * vmu;

mSig0i  = (eye(np) - mPhi.^2) * mSigi; % //Sig0^(-1)

%for(i=0 ; i<ns ; i++){
for i = 1:ns

    if i == 1

        mSigh   = inv( amX(i,:)'*amG2inv(:,:,i)*amX(i,:)...
                    + mSig0i + mSigi * mPhi.^2 );
        vbh     = mSigh...
                    * ( amX(i,:)'*amG2inv(:,:,i)*my(i,:)' + mSig0i * vmu...
                    + mSigi * mPhi * (mba(i+1,:)' - vpm1) );

    elseif i > 1 && i < ns

        mSigh   = inv( amX(i,:)'*amG2inv(:,:,i)* amX(i,:) + mSigp);
        vbh     = mSigh...
                    * ( amX(i,:)'*amG2inv(:,:,i)*my(i,:)' ...
                    + mSigi...
                    * ( mPhi * (vba_1 + mba(i+1,:)') + vpm2 ) );

    else

        mSigh   = inv( amX(i,:)'*amG2inv(:,:,i)* amX(i,:)  + mSigi );
        vbh     = mSigh...
                    * ( amX(i,:)'*amG2inv(:,:,i)*my(i,:)'... 
                    + mSigi * ( mPhi * vba_1 + vpm1 ) );

    end

    % Draw candidate (using the mean and variance from above. These are
    % based on the ith mba
    vbn = vbh + chol(mSigh)' * randn(np, 1); 
    vbo = mba(i,:)';

    mSighi  = inv(mSigh);
    dhn     = -0.5 * log(det(mSigh));
    dho     = dhn;

    % evaluate proposal state (vbn) and current state (vbo) against
    % mean and covaraince from proposal distribution (vbh and mSigh)
    %dhn += -0.5 * (vbn - vbh) ' mSighi * (vbn - vbh);
    %dho += -0.5 * (vbo - vbh) ' mSighi * (vbo - vbh);
    dhn = dhn + -0.5 * (vbn - vbh)'* mSighi * (vbn - vbh);
    dho = dho + -0.5 * (vbo - vbh)'* mSighi * (vbo - vbh);

    % Check proposed state (vbn) against current treshold
    %if(sumc(fabs(vbn) .< vd)){
    if any(abs(vbn) < vd)
          % some of the state variables are below the threshold

          % Only use thos who are above
          mxh = amX(i,:).* (abs(vbn) >= vd)';

          % draw new mean and covariance based on mxh and mba, where
          % mba is the current state
          if i == 1

              mSigh     = inv( mxh'*amG2inv(:,:,i)*mxh...
                              + mSig0i + mSigi * mPhi.^2 );
              vbh       = mSigh...
                            * ( mxh'*amG2inv(:,:,i)*my(i,:)' + mSig0i * vmu...
                            + mSigi * mPhi * (mba(i+1,:)' - vpm1) );

          elseif i > 1 && i < ns

              mSigh     = inv( mxh'*amG2inv(:,:,i)*mxh + mSigp );
              vbh       = mSigh...
                            * ( mxh'*amG2inv(:,:,i)* my(i,:)'...
                            + mSigi...
                            * ( mPhi * (vba_1 + mba(i+1,:)') + vpm2) );

          else

              mSigh     = inv( mxh'*amG2inv(:,:,i)*mxh + mSigi );
              vbh       = mSigh...
                            * ( mxh'*amG2inv(:,:,i)*my(i,:)'...
                            + mSigi * ( mPhi * vba_1 + vpm1 ) );

          end

          dln =	-0.5 * log(det(mSigh))...
                -0.5 * (vbn-vbh)'*inv(mSigh)*(vbn-vbh);

    else
        % none of the states are below the threshold, and dln is equal
        % to dhn
        dln = dhn;

    end

    % Check current state (vbo) against current treshold
    %if(sumc(fabs(vbo) .< vd)){
    if any(abs(vbo) < vd)

          mxh = amX(i,:) .* (abs(vbo) >= vd)';

          if i == 1

              mSigh     = inv( mxh'*amG2inv(:,:,i)*mxh...
                              + mSig0i + mSigi * mPhi.^2 );
              vbh       = mSigh...
                            * ( mxh'*amG2inv(:,:,i)* my(i,:)' + mSig0i * vmu...
                            + mSigi * mPhi * (mba(i+1,:)' - vpm1));

          elseif i > 1 && i < ns

              mSigh     = inv( mxh'*amG2inv(:,:,i)*mxh + mSigp );
              vbh       = mSigh...
                            * ( mxh'*amG2inv(:,:,i) * my(i,:)'...
                            + mSigi...
                            * ( mPhi * (vba_1 + mba(i+1,:)') + vpm2) );

          else

              mSigh     = inv( mxh'*amG2inv(:,:,i)*mxh  + mSigi );
              vbh       = mSigh...
                            * ( mxh'*amG2inv(:,:,i)* my(i,:)'...
                            + mSigi * ( mPhi * vba_1 + vpm1 ) );
          end

          dlo =	-0.5 * log(det(mSigh))...
                -0.5 * (vbo-vbh)'* inv(mSigh)*(vbo-vbh);

    else
        % none of the states are below the threshold, and dlo is equal
        % to dho            
        dlo = dho;

    end

    % calculate the acceptance probability
    dfrac = exp(dln - dhn - dlo + dho);

    if rand < dfrac
        mba(i,:)    = vbn';
        vba_1       = vbn;
        %ca++;
        ca          = ca + 1;
    else
        vba_1       = mba(i,:)';
    end

 % end of for loop   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /*
% **  fSampLT samples latent threshold
% **  conditional on time-varying parameter
% **
% **  [NOTE]
% **    y_t = X_t*b_t + e_t,    e_t ~ N(0, G2)
% */
% LTVAR::fSampLT(const my, const amX, const mp,
% 			   const amG2inv, const vd, const vdc)
% {
% 	decl ns, nk, vdn, vdo, mdn, mdo, dln, dlo;
% 	decl ve, ca, i, t;

function [vdo, ca] = drawBThreshold(my, amX, vd, mp, amG2inv, vdc)

% 	ns = rows(mp);
% 	nk = columns(mp);
% 	vdo = vd;
% 	ca = 0;

ns  = size(my,1);    	
nk  = size(amX,2);
% number of accepted draws
ca  = 0;
% current equals new
vdo = vd;


	
%for(i=0 ; i<nk ; i++){
for i = 1 : nk

    vdn     = vdo;
    vdn(i)  = vdc(i,i)*rand;        
      

%     mdn = fabs(mp) .< vdn' .? 0 .: mp;
%     mdo = fabs(mp) .< vdo' .? 0 .: mp;
    
    idx             = abs(mp) < repmat(vdn',[ns 1]);             
    vb              = vec(mp);
    vb(vec(idx))    = 0;
    mdn             = reshape(vb,[ns nk]);    
    
    idx             = abs(mp) < repmat(vdo',[ns 1]);             
    vb              = vec(mp);
    vb(vec(idx))    = 0;
    mdo             = reshape(vb,[ns nk]);        
    
    [dln,dlo] = deal(0);        
    for t = 1 : ns
%         ve = my[t][] - mdn[t][] * amX[t]';
%         dln += ve * amG2inv[t] * ve';
%         ve = my[t][] - mdo[t][] * amX[t]';
%         dlo += ve * amG2inv[t] * ve';

         ve = my(t) - mdn(t,:)*amX(t,:)';
         dln = dln + ve * amG2inv(:,:,t)*ve';
         
         ve = my(t) - mdo(t,:)*amX(t,:)';
         dlo = dlo + ve * amG2inv(:,:,t) * ve';

    end

    if rand < exp( -0.5 * (dln - dlo))
        vdo = vdn;
        ca  = ca + 1;
    end
    
end
