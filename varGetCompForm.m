function [BETAC, ALFAC] = varGetCompForm(beta, alfa, nlag, nvar)
% PURPOSE: Represent VAR(p) as VAR(1), eg comp form
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   alfa = vector. (n x 1), where n is the number of variables. Constants
%   form regression
%
%   beta = matrix. (n x (m-1)), where n is the number of variables, and m-1
%   is the number of variables*number of lags. 
%
%   nlag = number of lags in VAR 
%
%   nvar = number of variables in VAR
%
% Output:
%   betaC = matrix. (n+(n*(nlag-1))) x (m-1)). 
%
%   alfaC = vector. (n+(n*(nlag-1))) x 1). 
%
% Usage:
%   [BETAC,ALFAC]=varGetCompForm(beta,alfa,nlag,nvar) 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

[~,m1,d]=size(beta);
BETAC=nan(m1,m1,d);
if ~isempty(alfa)
    ALFAC=nan(m1,1,d);
else
    ALFAC=[];
end

for i=1:d
    if nlag>1
        BETAC(:,:,i)=[beta(:,:,i);
            eye((nvar)*(nlag-1)) zeros((nvar)*(nlag-1),nvar)];    
        if ~isempty(alfa)
            ALFAC(:,:,i)=[alfa(:,i);    
                zeros((nvar)*(nlag-1),1)];        
        end;
    else
        BETAC(:,:,i)=beta(:,:,i);
        if ~isempty(alfa)
            ALFAC(:,:,i)=alfa(:,i);    
        end;            
    end;                
end;
