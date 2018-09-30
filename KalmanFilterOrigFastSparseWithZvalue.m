function [state, varCov] = KalmanFilterOrigFastSparseWithZvalue(y, Z, D, H, T, C, R, Q, a0, p0) % #codegen
% KalmanFilterOrigFastSparseWithZvalue 
%
%   Faster version of Kalman Filter using sparse matrices as input, and
%   Zvalue object for z
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
num_observations    = size(y, 2);
num_states          = size(a0,1);

% Pre-allocate space for state and covariance of state
state   = NaN(num_states, 1, num_observations);
varCov  = NaN(num_states, num_states, num_observations);

a = a0;
P = p0;
I = eye(size(P, 1));
% 3. Subsequent observations
for t = 1 : 1 : num_observations
    
    Zt      = sparse(Z(:, :, t));     
    Tt      = sparse(T(:, :, t));
    Qt      = sparse(Q(:, :, t));
    Rt      = sparse(R(:, :, t));
    
    apt = Tt*a + C;
    Ppt = Tt*P*Tt' + Rt*Qt*Rt';     
    
    mask    = ~isnan(y(:, t));
    W       = diag(mask);
    Wt      = sparse(W(mask, :));
        
    if ~isempty(Wt)
        Ft = Wt * (Zt * Ppt * Zt' + H) * Wt';            
        
        y(~mask,t)  = 0; 
        vt          = sparse(Wt * ( y(:,t) - Zt * apt - D));
                        
        K = (Ppt * Zt' * Wt') / (Ft);
            
        a = apt + K * vt;    
        P = (I - K * Wt * Zt) *  Ppt;                        
        
    else 
        a = apt;
        P = Ppt;
    end
    
    state(:, 1, t)  = a;
    varCov(:, :, t) = P;
    
end




