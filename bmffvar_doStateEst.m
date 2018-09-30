function state = bmffvar_doStateEst(y, z, d, h, t, c, r, q, a0, p0, p, a0z,...
                                        numDynFac, numStates, indexes, initial)                                                    

nFreq   = ~isempty(indexes.qIndex)*2 + ~isempty(indexes.mIndex)*2 + ~isempty(indexes.wIndex)*2;

[T, N]  = size(y);
% make zadj to be used below
zadj    = cat(2, z, zeros(N, nFreq, T));                            

% Adjust y with AR(1) residuals
y    = y - repmat(d', [T 1]);
yadj = getWhiteningY(y, p);                  

% make z one observation longer since we loose one observation when making
% AR(1) correction
z   = cat(3, a0z, z); 

zi  = getWhiteningZTVP(z, p, numDynFac, nFreq + 1);                             
zadj(indexes.z_dIndex_r, indexes.z_dIndex_c : indexes.z_dIndex_c + 1, :)  = zi(indexes.z_dIndex_r, [indexes.z_dIndex_c indexes.z_dIndex_c*2], :);                                
zadj([indexes.z_qIndex_r(:);indexes.z_mIndex_r(:);indexes.z_wIndex_r(:)], end - nFreq + 1 : end, :)   = zi([indexes.z_qIndex_r(:);indexes.z_mIndex_r(:);indexes.z_wIndex_r(:)], indexes.z_dIndex_c + 1 : indexes.z_dIndex_c + nFreq, :);    


tadj = cat(1, t, zeros(nFreq, nFreq + numStates, T));
tadj = cat(2, tadj, zeros(nFreq*2 + numStates, nFreq, T));
tadj(end - nFreq + 1 : end, 1 : nFreq, :) = repmat(eye(nFreq), [1 1 T]);

radj = cat(1, r, zeros(nFreq, 1, T));
cadj = sparse(cat(1, c, zeros(nFreq, 1)));

d       = sparse(zeros(N, 1));
h       = sparse(h);
a0adj   = sparse(cat(1, a0, zeros(nFreq, 1)));
p0adj   = sparse(blkdiag(p0, eye(nFreq).*100));
    
%% Do filter
[state1, varCov1] = KalmanFilterOrigFastSparseWithZvalue(yadj', zadj,...
                                                 d, h, tadj, cadj, radj, q,...
                                                 a0adj, p0adj);
state1  = state1(1 : nFreq + numStates, 1, :);
varCov1 = varCov1(1 : nFreq + numStates, 1 : nFreq + numStates, :);
                                            
qq      = cat(3, initial.Qirep, full(q));
rr      = cat(3, initial.Rirep, full(r));
try
    QQ = mtimesx(mtimesx(rr, qq), rr, 'T');
catch
    QQ = zeros(size(r, 1), size(r, 1), T + numStates);
     for i = 1 : T + numStates
         QQ(:, :, i) = rr(:, :, i)*qq(:, :, i)* rr(:, :, i)';
     end  
end

a0vec   = repmat(a0, [1 1 numStates]);
p0vec   = repmat(p0, [1 1 numStates]);
state2  = cat(3, a0vec, state1(:, 1, :));
varCov2 = cat(3, p0vec, varCov1);
t       = cat(3, initial.Tirep0, t);

J                       = 1 + nFreq;
obsIdx                  = false;
state                   = kalmanFilterBackwardStateRecursionTVP(...
                            state2, varCov2, J, t, QQ(1 : 1 + nFreq, 1 :  1 + nFreq, :), obsIdx);   
state(:, 1 : nFreq)     = permute(state2(1 : nFreq, 1, :), [3 1 2]);                  
                  
