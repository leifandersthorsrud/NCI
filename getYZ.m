function [y,z,k]=getYZ(y,nlag,constant,betaRestrictions)
% PURPOSE: Put VAR on a form where vec(Y)= B*Z+e, such that different
% equations can contain different regressors
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

nvary       = size(y,2);          
[esty, estx] = varOrderData(y, nvary, nlag, constant);    

z = [];
% K is empty storage for number of variables in each eqaution
k = [];
for i = 1:nvary
    z = blkdiag(z, estx(:,betaRestrictions(i,:) == 1));
    %Z=kron(eye(neqs),estX);    
    k = cat(1, k, sum(betaRestrictions(i,:) == 1));
end
y = esty(:);