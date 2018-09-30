function BB=kalmanFilterBackwardStateRecursionTVP(B,P,J,F,Q,obsIdx,C)
% PURPOSE: Do backward recursion of state in Kalman Filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%
%   B = array, (ns x 1 x t), where ns is number of states, and t number of
%   observations (state output from KF iterations)
%
%   P = array, (ns x ns x t), covariance of state (varCov output from KF 
%   iterations)
%   
%   J = vector (n x 1), where n is the number of states to consider in the
%   backward recursions, i.e. when the Q matrix (see below) is singular
%   (e.g. when lags are included in the state) (n<=ns)
%
%   F = matrix (ns x ns x t). Coeff matrix for transition equation in KF
%
%   Q = matrix (n x n x t). Covariance of transition equaition residuals 
%
%   obsIdx = vector with logicals (n x 1). Ones for observable variables in 
%   the state equation (these should not be "updated", they are known), 
%   zeros for non-observables. 
%
%   C = vector (ns x 1 x t) with constants for the transition equation. This
%   input can be dropped, in which case it is assumed to be zeros.
%
% Output:
% 
%   BB = matrix (n x t) with "updated" states. 
%
% Usage:
%
%   BB=kalmanFilterBackwardStateRecursion(B,P,J,F,Q,obsIdx)
%
% Note: No error checking of input to ensure speed!
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
         

    T = size(B, 3);          

    if nargin == 6
        C = zeros(numel(J), 1, T);
    else
        C = C(J, 1, T);
    end;
      
    BB = nan(numel(J), T);                    
    Jl = J( obsIdx == 0 );            
    Jo = J( obsIdx == 1 );               
    
    try
        BB(Jl, T) = B(Jl, 1, T) + chol(P(Jl, Jl, T))'*randn(numel(Jl), 1);
        BB(Jo, T) = B(Jo, 1, T);    
    catch
        BB(J, T) = B(J, 1, T);
    end
    
    for t = T-1:-1:1  
        Pt 	= P(:,:,t);
        Ft 	= F(J,:,t);
        Bt  = B(:,1,t);

        m               = ( (Pt*Ft')/(Ft*Pt*Ft' + Q(J,J,t)) );                        
        bb              = Bt + m*(BB(J,t+1) - C(:,1,t) - Ft*Bt);
        pp              = Pt - m*Ft*Pt;                                
        try
            BB(Jl, t)   = bb(Jl) + chol(pp(Jl, Jl))'*randn(numel(Jl), 1);
            BB(Jo, t)   = bb(Jo);               
        catch Me
            BB(J, t)   = Bt(J);  
                m           = ( (Pt(Jl, Jl)*F(Jl,Jl,t)')*pinv(F(Jl,:,t)*Pt(Jl,Jl)*F(Jl,:,t)' + Q(Jl,Jl,t)) );                        
                bb          = Bt(Jl) + m*(BB(Jl,t+1) - C(Jl,1,t) - F(Jl,:,t)*Bt);
                pp          = Pt(Jl, Jl) - m*F(Jl,Jl,t)*Pt(Jl, Jl);                                                        
                BB(Jl, t)   = bb + chol(pp)'*randn(numel(Jl), 1);                
        end
                         
    end
    
    BB = BB';
        
end