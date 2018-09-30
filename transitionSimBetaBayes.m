function [T,C] = transitionSimBetaBayes(nobs, nlag, constant, ...
    ks, ys, zs, hs, Vprior, Tprior, Tres)
% transitionSimBayes   -  Draw parameters from  state equation using
%                          Bayesian methods
%
% Input:
%
% Output:
%
% Usage:
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

maxNumberOfRhs  = max(ks);
nvary           = numel(ks);

stabilityCondition = 1;
while stabilityCondition == 1
    
    % Compute some usefull parameters    
    xhx         = zs'*hs*zs;
    xhy         = zs'*hs*ys;
    capv0inv    = inv(Vprior);
    
    %draw from beta conditional on H
    capv1inv    = capv0inv + xhx;
    capv1       = inv(capv1inv);
    b1          = capv1*(capv0inv*Tprior + xhy);
    beta        = b1 + norm_rnd(capv1);    
        
    % Transform ALPHA into form that can be used in other utils below and
    % stored as output
    betan = zeros(nvary, maxNumberOfRhs);
    cnt = 0;
    for i = 1:nvary
        betan(i, Tres(i,:) == 1) = beta(1 + cnt:ks(i) + cnt);
        cnt = cnt + ks(i);
    end;
    
    % get comp form to check stability
    if constant
        [T,C] = varGetCompForm(betan(:,2:end,:),betan(:,1,:), nlag, nvary);
    else
        [T,C] = varGetCompForm(betan, [], nlag, nvary);
    end;
    
    if any(abs(eig(T)) >= 1)
        stabilityCondition  = 1;%1
    else
        % To brake out of while statement if stability
        % ensured
        stabilityCondition  = 0;        
    end
end