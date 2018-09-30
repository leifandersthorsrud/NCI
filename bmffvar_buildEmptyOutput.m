function [TS, QS, HS, AS, ZS, WS, PS, ddS, RS, DS, BS] = ...
            bmffvar_buildEmptyOutput(numStatesL, nObs, N, nsave,...
            numDynFac, indexes)

nq = numel(indexes.qIndex);        
nm = numel(indexes.mIndex);        
nw = numel(indexes.wIndex);        
nd = numel(indexes.dIndex);

nvaraL = indexes.nvaraL;

TS = ndSparse.build([numStatesL, numStatesL, nObs, nsave]);
QS = ndSparse.build([1, 1, nObs, nsave]);
HS = ndSparse.build([N, N, 1, nsave]);
AS = ndSparse.build([numStatesL, 1, nObs, nsave]);
ZS = ndSparse.build([N, numStatesL, nObs, nsave]);
PS = ndSparse.build([N, 1, 1, nsave]);
ddS = ndSparse.build([N, 1 + numDynFac, 1, nsave]);
RS = ndSparse.build([numStatesL, 1, nObs, nsave]);
DS = ndSparse.build([N, 1, 1, nsave]);
BS = ndSparse.build([1, 1, 1, nsave]);

N = nq + nm + nw + nd*(numDynFac + 1);
WS = ndSparse.build([N, N, 1, nsave]);


