function [sigmadraw, bdraw, s] = drawSigmaTVP(y, sigmadraw, bdraw,...
                                                h0_mean, h0_cov, vprior, Sprior)
% PURPOSE: Draw stochastic volatility states and the variance of the
% variance. Useing Kim et. al (1998) procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Model: 
% 
% y=Sigma * u  where u ~ N(0,1) and Sigma is a diagonal matrix with std's
% on the diagonal. By squaring terms and taking logs, we get: 
%
% y*=2h(t) + e(t)    where e=log(u^2), h(t)=log(Sigma) and y*=log(y^2)        
% h(t)=h(t-1) + b(t) where the covariance of b(t) is B 
%
% Note 1) that the covariance of e(t) is an identity matrix. This is why we
% can use the Kim et. al. procedure below.  
%
% Note 2) that the y vector should not be on the y* form, but only y. The
% transformation is done within the code, adding an offsetting constant. 
%
% Input:
%
% y = (t x n) matrix with observables (transfomation to logs and adding a
% offsetting constant is done inside the function
%
% sigmadraw = (t x n) matrix with sigmadraw from earlier recursion
% of the function. Used to update draw the s index. Important, these should
% be exp(h(t)), i.e., standard deviation. Inside the code we transform them
% to log variances. 
%
% bdraw = (n x n) matrix with previous draws for covariance of the
% variances. 
%
% h0_mean = (n x 1) vector with a0 values. These should in most cases
% be 0's. 
%
% h0_cov = (n x n) matrix of p0 values. Typically an identity matrix
% scales up by a constant (to make it more diffuse)
%
% vprior = scalar. Degrees of freedom parameter. Used for Wishart
% dirstibution. Must exceed the dimension of the W to make the prior
% proper. 
%
% Sprior = (n x n) matrix. Prior for the (log) variance of the variance. 
%
% Output: 
%
% s = (t x n) matrix with indicators from earlier recursion
% of the function. Used to mix the mixture of draws. Must have elements
% with values in the range 1,2,...7. Updated/sampled from input based on
% simulated values within the function. 
%
% bdraw = new draw of the covariance of the variance. See above
%
% sigmadraw = (t x n) matrix. Draw of state vector, but scaled and 
% transformed such that they are std's, i.e., exp(h(t)). 
%
% Note 3): 
%
% The function only works for the model as specified above. I.e., random
% walk assumption for the volatilities. The y input matrix should be
% transformed appropriately before used as input. 
%
% Usage: 
%
% [sigmadraw,bdraw,s]=drawSigmaTVP(y,sigmadraw,bdraw,h0_mean,h0_cov,vprior,Sprior)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Copyright (C) 2014  PROFOR Team
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get some sizes
[Tt, n] = size(y);

% Parameters of the 7 component mixture approximation to a log(chi^2)
% density:
q_s     = [0.00730; 0.10556; 0.00002; 0.04395; 0.34001; 0.24566; 0.25750];     % probabilities
m_s     = [-10.12999; -3.97281; -8.56686; 2.77786; 0.61942; 1.79518; -1.08819];% means
u2_s    = [5.79596; 2.61369; 5.17950; 0.16735; 0.64009; 0.34023; 1.26261];    % variances
num_m_s = numel(m_s);

% Off-setting constant
c       = 0.001;
% Construct yss= log(y^2 + c)
yss     = log(y.^2 + c);

% First draw s (chi square approximation mixture component) conditional on 
% everything else. This is used to update at the next step the log-volatilities 
% Sigtdraw
Sigtdraw    = log(sigmadraw.^2);
s           = nan(Tt, n);
for jj = 1:n
    for i = 1:Tt
        prw = nan(num_m_s, 1);
        for j = 1:num_m_s
            % NOTE: Using yss
            %temp1=(1/sqrt(2*pi*u2_s(j)))*exp(-.5*(((yss(i,jj) - 2*Sigtdraw(i,jj) - m_s(j) + 1.2704)^2)/u2_s(j)));
            temp1       = ( 1/sqrt(2*pi*u2_s(j))) * exp( -.5*((( yss(i,jj) - Sigtdraw(i,jj) - m_s(j) + 1.2704)^2 )/u2_s(j)) );
            prw(j,1)    = q_s(j)*temp1;
        end
        prw     = prw./sum(prw);
        cprw    = cumsum(prw);
        trand   = rand(1, 1);
        if trand < cprw(1, 1); 
            imix = 1;
        elseif trand < cprw(2, 1)
            imix = 2;
        elseif trand < cprw(3, 1)
            imix = 3;
        elseif trand < cprw(4, 1)
            imix = 4;
        elseif trand < cprw(5, 1)
            imix = 5;
        elseif trand < cprw(6, 1)
            imix = 6;
        else
            imix = 7;
        end
        s(i,jj) = imix;  % this is a draw of the mixture component index
    end
end

% In order to draw the log-volatilies, substract the mean and variance
% of the 7-component mixture of Normal approximation to the measurement
% error covariance
H       = zeros(n, n, Tt);
yss1    = zeros(Tt, n);
for i = 1:Tt
    for j = 1:n
        imix        = s(i, j);
        H(j, j, i)  = u2_s(imix);                
        yss1(i, j)  = yss(i, j) - m_s(imix) + 1.2704;
    end
end

% Do KF and state draws
% NOTE: Using yss1
Sigtdraw = kalmanFilterStep(n, Tt, H, bdraw, h0_mean, h0_cov, yss1);
% Sigtdraw is a draw of the diagonal log-volatilies, which will give us
% SIGMA(t), or h(t) in the notation used above

% Draws in Sigtdraw are in logarithmic scale (log-volatilies). Create 
% original standard deviations of the VAR covariance matrix
%sigmadraw=exp(Sigtdraw);
sigmadraw = exp(0.5.*Sigtdraw);

%Draw W, the covariance of SIGMA(t) (from iWishart)
emat    = Sigtdraw(2:end,:) - Sigtdraw(1:end-1,:);
v       = Tt + vprior;
H       = inv(vprior*Sprior + emat'*emat);    
bdraw   = inv(wish_rnd(H, v));        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = kalmanFilterStep(n, Tt, Ht, bdraw, a0, p0, yss1)

z = repmat(eye(n), [1 1 Tt]);  % Zs
d = zeros(n, 1, Tt);
h = Ht;  % vart
t = repmat(eye(n), [1 1 Tt]);
c = zeros(n, 1, Tt);
r = repmat(eye(n), [1 1 Tt]); 
q = repmat(bdraw, [1 1 Tt]);  % bdraw

% Do KF filter
[state, varCov] = KalmanFilterOrigFastNoNanTVP(yss1', z, d, h, t,...
                                                                c, r, q, a0, p0);       

% Do backwards recursion to generate state                 
J       =(1:n);
obsIdx  = false(n, 1);
state   = kalmanFilterBackwardStateRecursion(state, varCov,...
                                                J, t(:,:,1), q(:,:,1), obsIdx);                               
