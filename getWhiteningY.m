function ystar = getWhiteningY(y, p)
% getWhiteningY - clean y for AR components in residuals
%
% Input: 
%
%   y - matrix (t x n), where t is number of observations and n number of
%   series
%
%   p - matrix (n x m), where m is number of lags in AR processes
%
% Output: 
%
%   ystar - matrix (t x n)
%
% Usage: 
%
%   ystar = getWhiteningY(y, p)
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

[T, n]  = size(y);
nlag    = size(p, 2);
ystar   = y;
ynlag   = latMlag(y, nlag, 0);               

nanIdx        = isnan(ynlag);
ynlag(nanIdx) = 0;

cnt = 0;
for i = 1 : nlag
    ystar = ystar - ynlag(:, 1 + cnt:n + cnt).*repmat(p(:, i)',[T 1]);
    cnt = cnt + n;
end