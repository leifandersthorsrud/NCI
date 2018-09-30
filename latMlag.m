function xlag = latMlag(x, p, init)
% latMlag - Generate lags of matrix (vector), following common VAR standards, eg: 
%           x1(t-1) x2(t-1) ... xN(t-1)... x1(t-p)... xN(t-p). 
%   NOTE: This is NOT as done in the Lesage package, which uses x1(t-1) ... 
%   x1(t-p) ... xN(t-1)...xN(t-1)..xN(t-p)
%
% Input:
%   x = matrix or vector, (t x n) where t is number of obs, and n is number
%   of variables.
%
%   p = integer. Number of lags. 
%
%   init = (optional) scalar value to feed initial missing values
%   (default = 0)
%
% Output:
%   xlag = matrix or vector with lags of x. The output has the same size as
%   the input matrix or vector. Missing values at the beginning (due to
%   taking lags) have been replaced by the init argument. 
%
% Usage:
%   xlag = latMlag(x,p,init)
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

if nargin ==1 
    p = 1; % default value
    init = 0;
elseif nargin == 2
    init = 0;
end

if nargin > 3
    error('mlag:err','Wrong # of input arguments');
end

[nobs, nvar] = size(x);

xlag = ones(nobs,nvar*p)*init;

for ii=1:p
    xlag(1+ii:nobs,(nvar*(ii-1)+1):nvar*ii)=x(1:nobs-ii,:);
end
