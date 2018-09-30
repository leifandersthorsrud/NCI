function [estY,estX,vname,dates,xtilde]=varOrderData(data,nvar,nlag,constant)
% PURPOSE: Get right hand side variables and lefthand side variables.
% removes any rows with nans and adjusts any dates acordingly 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   data = matrix. (t x n), where t is the number of observations in the
%   sample, and n is the number of variables to be used in the VAR
%   estimation
%
%   nlag = number of lags in VAR 
%
%   nvar = number of variables in VAR
%
%   varargin(optional)
%   'vnames' = string, followed by a cell (1 x n) with the names of the
%   variables in the data matrix, eg {'gdp','inf'}, in the same order as
%   the data in the data matrix
%   
%   'dates' = string, followed by a vector of CS dates (t x n). 
%
%   'constant' = string, followed by a logical. If true constant included
%   in regression. Default=true, else false. 
%
% Output:
%   estY = matrix. (t x n), where t is the number of observations, and n is
%   the number of variables in the VAR
%
%   estX = matrix. (t x m), where t is the number of observations, and m is
%   the number of variables*number of lags + 1 (constant). Note that the
%   ordering of the variables should be according to lags (and not
%   variables), eg y(t-1), y(t-2) etc. 
%
%   vname = cell string. (1,1+n*nlags), with variable names. Matches the
%   ordering of the estX variables. If vnames provided in varargin, these
%   will be used, if not defualts=coeff1..coeffn..
%
%   dates = vector with CS dates. (tt x1), where tt is the length of the
%   estY matrix. If dates provided in varargin these will be truncated and
%   used, if not this will be a be vector of nans
%
%   xtilde = vector, (1 x nvar*nlag + constant). Vector that can be used to 
%   generate forecasts. If constant in model, constant=1, else 0. This will
%   be ordered as ordering in estY and with number of lags provided. See
%   also estX

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

[t,c]=size(data);
vname=makeVnames(c);

options.dates       = (1:t)';
options.vname       = vname;
options.constant    = constant;

data=[options.dates data];
% remove nans

data(any(isnan(data),2),:)=[];    
% pull out different parts
dates=data(:,1);
y=data(:,2:end);
const=ones(size(data,1),1);
% set up estimation eq
xmat=latMlag(y,nlag,NaN);

if options.constant
    estData=[dates y const xmat];           
else    
    estData=[dates y xmat];
end;

% remove nans
estData(any(isnan(estData),2),:)=[];
% make rest of output
estY=estData(:,2:nvar+1);
estX=estData(:,nvar+2:end);

% variable names   
vname=repmat(options.vname(1:end),[1 nlag]);
if options.constant    
    vname=[{'constant'}, vname];    
end;    
dates=estData(:,1);

tmp=estY(end:-1:end-nlag+1,:)';
if options.constant        
    xtilde=[1 tmp(:)'];            
else
    xtilde=tmp(:)';
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vname=makeVnames(c)
vname=cell(1,c);
for n=1:c
    vname{1,n}=['coeff' int2str(n)];
end;
