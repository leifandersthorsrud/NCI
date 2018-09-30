%% This file replicates the simulation experiment done in online Appendix C.1. 
% of the article:
% Thorsrud, Leif Anders, 2018. "Words are the new numbers. A newsy coincident index of the
% business cycle," Journal of Business & Economic Statistics, forthcoming
%
% If you use these codes, please use the provided reference 
% (Note, the codes are not written for efficiency nor elegance. The codes 
% are just meant for simple replication)
%

%% 1) Load data

% The data.mat contains two objects; 'dates' and 'y'. The 'y' matrix
% contains the simulated data. The first column is the quarterly values,
% the second column is the monthly values, the third column is the weekly
% values, and the remaining column are the daily values. Note that the
% quarterly, monthly, and weekly values contain nans. Observations are only
% included at the assumed observation interval, i.e., the last day of the
% period. The dates vector has equal length as the y matrix
data = load('data.mat');
% Extract
y = data.y;

%% 2) Load state-space matrices with deterministic elements filled in (see 
%% article for details + helper objects to keep track of the different matrix 
%% elements

% These elements will need to be changed if you want to estimate a model
% with different lag specification, number of daily, weekly, monthly, and
% quarterly variables etc. Thus, this script is just a simple replication
% of the simulation experiment done in Appendix C.1. Also, in the loaded SS
% matrices, starting values for the MCMCM are filled in. These can 
% be changed.
SSmatrices      = load('SSmatrices.mat');
% Extract
Cirep   = SSmatrices.Cirep;
Qirep   = SSmatrices.Qirep;
Tirep   = SSmatrices.Tirep;
Tirep0  = SSmatrices.Tirep0;
Rirep   = SSmatrices.Rirep;
Direp   = SSmatrices.Direp;
Hirep   = SSmatrices.Hirep;
Zirep   = SSmatrices.Zirep;
Zirepu  = SSmatrices.Zirepu;
Pirep   = SSmatrices.Pirep;
dtIrep  = SSmatrices.dtIrep;
Wirep   = SSmatrices.Wirep;
Birep   = SSmatrices.Birep;

a0A     = SSmatrices.a0A;
p0A     = SSmatrices.p0A;
a0Z     = SSmatrices.a0Z;
initial = SSmatrices.initial;

% These elements will need to be changed if you want to estimate a model
% with different lag specification, number of daily, weekly, monthly, and
% quarterly variables etc. Thus, this script is just a simple replication
% of the simulation experiment done in Appendix C.1 
helperObjects   = load('helperObjects.mat');
% Extract
indexes             = helperObjects.indexes;
nlag                = helperObjects.nlag;
numStates           = nlag;
numStatesL          = 2*3 + numStates;
numDynFac           = 0;
nObs                = helperObjects.nObs;
N                   = helperObjects.N;
observableIdx_dq    = helperObjects.observableIdx_dq;
observableIdx_dm    = helperObjects.observableIdx_dm;
observableIdx_dw    = helperObjects.observableIdx_dw;

TRes1               = ones(1, nlag);

%% 3) Load the priors (change them if you want)

priorSettings = load('priors.mat');
priorSettings = priorSettings.priorSettingsStruct;

%% 4) Do the Gibbs MCMCM

% First element is number of burninn, second element is number of draws to
% save (high memory consumption...)
nSaveDraws = [1000 500];
nburn      = nSaveDraws(1);
% Number of steps to keep from the MCMC chain
nstep      = 2;

%% 4.0) Empty output

[~, QS, ~, AS, ZS, ~, ~, dtS, ~, ~, ~] = ...
            bmffvar_buildEmptyOutput(numStatesL, nObs, N,...
            nSaveDraws(2), numDynFac, indexes);    
                        
for di = 1 : numel(nSaveDraws)
    
    ntoti   = nSaveDraws(di)*nstep;        
    
    cnt = 1;    
    for irep = 1 : ntoti                                    
    
        % Disp irep                
        disp(irep)
        
        %% 4.1) State    
        Airep       = bmffvar_doStateEst(y, Zirepu, Direp, Hirep,...
                                        Tirep, Cirep, Rirep, Qirep, a0A, p0A,...
                                        Pirep, a0Z, numDynFac, numStates,...
                                        indexes, initial);          
        AirepLarge  = [Airep latMlag(Airep(:, end), nlag - 1)];                                                           

        %% 4.2) Observation equation
        %% Do: Pq, Hq, Zq Zuq Wq dq
        if ~isempty(indexes.qIndex)
            AirepLargeL = AirepLarge(nlag + 1 : end, indexes.z_qIndex_c);            
            AirepLargeL = [0;
                           AirepLargeL(observableIdx_dq, :)];
            AirepLargeL = [AirepLargeL latMlag(AirepLargeL, nlag - 1)];                        
            [p, h, z, zu, w, d] = bmffvar_doObsEst(...
                            y(observableIdx_dq, indexes.z_qIndex_r) - repmat(Direp(indexes.z_qIndex_r)', [sum(observableIdx_dq) 1]),...
                            AirepLargeL,... % a
                            Pirep(indexes.z_qIndex_r),...
                            Zirep(indexes.z_qIndex_r, indexes.z_qIndex_c, observableIdx_dq),...
                            dtIrep(indexes.z_qIndex_r, :),...
                            Hirep(indexes.z_qIndex_r, indexes.z_qIndex_r),...
                            Wirep(indexes.qIndex, indexes.qIndex),... 
                            priorSettings,...
                            numDynFac,...
                            indexes.z_qIndex_r,...
                            true,...true
                            'q'...
                            );
            Zirep(indexes.z_qIndex_r, indexes.z_qIndex_c, :)    = bmffvar_z2day(z, observableIdx_dq);
            Zirepu(indexes.z_qIndex_r, indexes.z_qIndex_c, :)   = bmffvar_z2day(zu, observableIdx_dq);
            Wirep(indexes.qIndex, indexes.qIndex)               = w;     
            Pirep(indexes.z_qIndex_r)                           = p;
            Hirep(indexes.z_qIndex_r, indexes.z_qIndex_r)       = h;
            dtIrep(indexes.z_qIndex_r, :)                       = d;                        
        end
            
        %% Do: Pm, Hm, Zm Zum Wm dm
        if ~isempty(indexes.mIndex)
            AirepLargeL = AirepLarge(nlag + 1 : end, indexes.z_mIndex_c);            
            AirepLargeL = [0;
                           AirepLargeL(observableIdx_dm, :)];
            AirepLargeL = [AirepLargeL latMlag(AirepLargeL, nlag - 1)];
            [p, h, z, zu, w, d] = bmffvar_doObsEst(...
                            y(observableIdx_dm, indexes.z_mIndex_r) - repmat(Direp(indexes.z_mIndex_r)', [sum(observableIdx_dm) 1]),...
                            AirepLargeL,... % a
                            Pirep(indexes.z_mIndex_r),...
                            Zirep(indexes.z_mIndex_r, indexes.z_mIndex_c, observableIdx_dm),...
                            dtIrep(indexes.z_mIndex_r, :),...
                            Hirep(indexes.z_mIndex_r, indexes.z_mIndex_r),...
                            Wirep(indexes.mIndex, indexes.mIndex),... 
                            priorSettings,...
                            numDynFac,...
                            indexes.z_mIndex_r,...
                            true,...
                            'm'...
                            );
            Zirep(indexes.z_mIndex_r, indexes.z_mIndex_c, :)     = bmffvar_z2day(z, observableIdx_dm);
            Zirepu(indexes.z_mIndex_r, indexes.z_mIndex_c, :)    = bmffvar_z2day(zu, observableIdx_dm);
            Wirep(indexes.mIndex, indexes.mIndex)                = w;     
            Pirep(indexes.z_mIndex_r)                            = p;
            Hirep(indexes.z_mIndex_r, indexes.z_mIndex_r)        = h;
            dtIrep(indexes.z_mIndex_r, :)                        = d;                                             
        end

        %% Do: Pw, Hw, Zw Zuw Ww dw
        if ~isempty(indexes.wIndex)
            AirepLargeL = AirepLarge(nlag + 1 : end, indexes.z_wIndex_c);            
            AirepLargeL = [0;
                           AirepLargeL(observableIdx_dw, :)];
            AirepLargeL = [AirepLargeL latMlag(AirepLargeL, nlag - 1)];
            [p, h, z, zu, w, d] = bmffvar_doObsEst(...
                            y(observableIdx_dw, indexes.z_wIndex_r) - repmat(Direp(indexes.z_wIndex_r)', [sum(observableIdx_dw) 1]),...
                            AirepLargeL,... % a
                            Pirep(indexes.z_wIndex_r),...
                            Zirep(indexes.z_wIndex_r, indexes.z_wIndex_c, observableIdx_dw),...
                            dtIrep(indexes.z_wIndex_r, :),...
                            Hirep(indexes.z_wIndex_r, indexes.z_wIndex_r),...
                            Wirep(indexes.wIndex, indexes.wIndex),... 
                            priorSettings,...
                            numDynFac,...
                            indexes.z_wIndex_r,...
                            true,...
                            'w'...
                            );
            Zirep(indexes.z_wIndex_r, indexes.z_wIndex_c, :)     = bmffvar_z2day(z, observableIdx_dw);
            Zirepu(indexes.z_wIndex_r, indexes.z_wIndex_c, :)    = bmffvar_z2day(zu, observableIdx_dw);
            Wirep(indexes.wIndex, indexes.wIndex)                = w;     
            Pirep(indexes.z_wIndex_r)                            = p;
            Hirep(indexes.z_wIndex_r, indexes.z_wIndex_r)        = h;
            dtIrep(indexes.z_wIndex_r, :)                        = d;                                             
        end
        
        %% Do: Pd, Hd, Zd Zud Wd dd
        AirepLargeL = Airep(:, end);            
        AirepLargeL = [AirepLargeL latMlag(AirepLargeL, nlag - 1)];
        AirepLargeL = AirepLargeL(nlag : end, :);
        [p, h, z, zu, w, d] = bmffvar_doObsEst(...
                        y(:, indexes.z_dIndex_r) - repmat(Direp(indexes.z_dIndex_r)', [size(y, 1) 1]),...
                        AirepLargeL,... % a
                        Pirep(indexes.z_dIndex_r),...
                        Zirep(indexes.z_dIndex_r, indexes.z_dIndex_c, :),...
                        dtIrep(indexes.z_dIndex_r, :),...
                        Hirep(indexes.z_dIndex_r, indexes.z_dIndex_r),...
                        Wirep(indexes.dIndex, indexes.dIndex),... 
                        priorSettings,...
                        numDynFac,...
                        indexes.z_dIndex_r,...
                        true,... true
                        'd'...
                        );
        Zirep(indexes.z_dIndex_r, indexes.z_dIndex_c, :)     = z;
        Zirepu(indexes.z_dIndex_r, indexes.z_dIndex_c, :)    = zu;
        Wirep(indexes.dIndex, indexes.dIndex)                = w;     
        Pirep(indexes.z_dIndex_r)                            = p;
        Hirep(indexes.z_dIndex_r, indexes.z_dIndex_r)        = h;
        dtIrep(indexes.z_dIndex_r, :)                        = d;      

        %% 4.3) Transition equation: T, Q
        AirepLargeL             = Airep(:, end);
        [QirepS, Tirep2, Birep] = bmffvar_doTransEst(nlag, TRes1, AirepLargeL, Qirep(end, end, :),... 
                                    Birep, priorSettings);                        
        Qirep(end, end, :)      = QirepS;        
                         
        Tirep(end - numStates + 1 : end, end - numStates + 1 : end, :)  = repmat(Tirep2, [1 1 nObs]);      
        Tirep0(end - numStates + 1 : end, end - numStates + 1 : end, :) = repmat(Tirep2, [1 1 nObs]);
                       
        if ~isempty(indexes.qIndex)            
            Tirep(indexes.t_qIndex_r, indexes.ta_Index_c, :) = Tirep0(indexes.t_qIndex_r, indexes.ta_Index_c, :) .* repmat(Tirep2(1, :), [2 1 nObs]);                            
        end            
        if ~isempty(indexes.mIndex)           
            Tirep(indexes.t_mIndex_r, indexes.ta_Index_c, :) = Tirep0(indexes.t_mIndex_r, indexes.ta_Index_c, :) .* repmat(Tirep2(1, :), [2 1 nObs]);                                        
        end                                                
        if ~isempty(indexes.wIndex)
            Tirep(indexes.t_wIndex_r, indexes.ta_Index_c, :) = Tirep0(indexes.t_wIndex_r, indexes.ta_Index_c, :) .* repmat(Tirep2(1, :), [2 1 nObs]);                                        
        end                                                        


        %% 4.4) Collect 
        
        % After burinn (taking into account that only every "step" is stored)
        if irep > nburn*nstep  && rem(irep, nstep) == 0                                             
                         
             QS(:, :, :, cnt)       = single(Qirep);                        
             AS(:, :, :, cnt)       = permute(single(AirepLarge(nlag + 1 : end, :)), [2 3 1]);                          
             ZS(:, :, :, cnt)       = single(Zirep);                          
             dtS(:, :, 1, cnt)      = single(dtIrep);             

             cnt = cnt + 1;
        end    
                        
    end
           
    % Set nburn to zero, to signal that we are done with the burninn
    % sampling
    nburn = 0;    
            
end

%% 5) Save the output

save('results.mat',...
    'QS', 'AS', 'ZS', 'dtS');

%% 6) See the file Graphs.m for graphing the output