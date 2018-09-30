function [state,varCov]=KalmanFilterOrigFastNoNanTVP(y,Z,d,H,T,c,R,Q,a0,p0)
% PURPOSE:	Apply the Kalman Filter to a time-varying state space
% formulation, specified as in Harvey (1989). Faster version of 
% KalmanFilterOrig, i.e no error checking, and
% no analytic derivation option for derivatives. See KalmanFilterOrig for 
% input desc. See especially last Note. 
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

% 1. Initialisation
%[num_observables,num_observations]=size(y);
num_observations    = size(y,2);
num_states          = size(a0,1);

% Pre-allocate space for state and covariance of state
state                = NaN(num_states,1,num_observations);
varCov               = NaN(num_states,num_states,num_observations);

a = a0;
P = p0;
I = eye(size(P, 1));

% 3. Subsequent observations
for t = 1 : 1 : num_observations                             
    
    Tt = T(:,:,t);
    ct = c(:,:,t);
    Rt = R(:,:,t);
    Qt = Q(:,:,t);
    Zt = Z(:,:,t);
    dt = d(:,:,t);
    Ht = H(:,:,t);
    
    apt = Tt*a + ct;
    Ppt = Tt*P*Tt' + Rt*Qt*Rt';         
    
    Ft  = Zt * Ppt * Zt' + Ht;   
    vt  = y(:, t) - Zt * apt - dt;    
    K   = (Ppt * Zt') / (Ft);    
            
    a   = apt + K * vt;    
    P   = (I - K * Zt) *  Ppt;                            
    
    state(:, 1, t)  = a;
    varCov(:, :, t) = P;    
   
end




