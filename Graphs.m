%% This file replicates the output reported in online Appendix C.1. 
% of the article:
% Thorsrud, Leif Anders, 2018. "Words are the new numbers. A newsy coincident index of the
% business cycle," Journal of Business & Economic Statistics, forthcoming
%
% Before running this script, run the file called Main.m
%

%% 1) Load result and true simulated data
% First, we load some simulated data, stored as Tsdata objects.
load('results.mat');
load('bmffvarData_simulation.mat');

%%

lengthToPlot = 5000;
quants       = [0.05 0.5 0.95];

%%

estA  = permute(quantile(full(AS), quants, 4), [3 4 1 2]);
estQ  = permute(quantile(full(QS), quants, 4), [3 4 1 2]);
nobs  = size(estA, 1);
trueA = S.Airep(end - nobs + 1 : end, end);
trueQ = S.Qirep(end - nobs + 1 : end);

estA  = estA(1 : lengthToPlot, :, :);
estQ  = estQ(1 : lengthToPlot, :, :);
trueA = trueA(1 : lengthToPlot, :);
trueQ = trueQ(1 : lengthToPlot, :);

% Correlations
corr(trueA, estA(:, 2, 7))
corr(trueQ, estQ(:, 2))


estAs7 = (estA(:, :, 7) -  repmat(mean(estA(:, 2, 7)), [size(estA, 1) 3]))./repmat(std(estA(:, 2, 7)), [size(estA, 1) 3]);
estQs  = (estQ - repmat(mean(estQ(:, 2)), [size(estQ, 1) 3])) + repmat(mean(trueQ), [size(estQ, 1) 3]);


fh = figure('name', 'sim_factor_true_factor');
p1 = plot(estAs7(:, 2), 'color', [0 0 0], 'linestyle', '-', 'marker', 'none', 'linewidth', 0.5);
hold on
p2 = plot(estAs7(:, [1 3]), 'color', [0 0 0], 'linestyle', '--', 'marker', 'none', 'linewidth', 0.5);
p3 = plot(zscore(trueA), 'color', [0.5 0.5 0.5], 'linestyle', '-', 'marker', 'none', 'linewidth', 0.5);
scrsz = get(fh, 'Position');
set(fh, 'Position', [scrsz(1) scrsz(2) scrsz(3)*2 scrsz(4)*2/3]);                
lh = legend([p1 p3], {'Posterior', 'True'}, 'box', 'off', 'fontsize', 12, 'edgecolor', [1 1 1], 'location', 'northeast');
box off
axis tight


fh = figure('name', 'sim_sw_true_sw');
p1 = plot(estQs(:, 2), 'color', [0 0 0], 'linestyle', '-', 'marker', 'none', 'linewidth', 0.5);
hold on
p2 = plot(estQs(:, [1 3]), 'color', [0 0 0], 'linestyle', '--', 'marker', 'none', 'linewidth', 0.5);
p3 = plot(trueQ, 'color', [0.5 0.5 0.5], 'linestyle', '-', 'marker', 'none', 'linewidth', 0.5);
scrsz = get(fh, 'Position');
set(fh, 'Position', [scrsz(1) scrsz(2) scrsz(3)*2 scrsz(4)*2/3]);                
lh = legend([p1 p3], {'Posterior', 'True'}, 'box', 'off', 'fontsize', 12, 'edgecolor', [1 1 1], 'location', 'northeast');
box off
axis tight
set(gca, 'ylim', [0 2])

%% quarterly, monthly and weekly loadings

ZquantileAll = quantile(full(ZS), quants, 4);
ZquantileAll = ZquantileAll(:, :, 1 : lengthToPlot, :);

Zq = permute(ZquantileAll(1, 1, :, : ), [3 4 1 2]);
Zm = permute(ZquantileAll(2, 3, :, :), [3 4 1 2]);
Zw = permute(ZquantileAll(3, 5, :, :), [3 4 1 2]);

figure
subplot(3,1,1)
plot(Zq, 'b')
hold on
tmp = squeeze(S.Zirep(1, 1, end - nobs + 1 : end));
plot(tmp(1 : lengthToPlot), 'r')
subplot(3,1,2)
plot(Zm, 'b')
hold on
tmp = squeeze(S.Zirep(2, 3, end - nobs + 1 : end));
plot(tmp(1 : lengthToPlot), 'r')
subplot(3,1,3)
plot(Zw, 'b')
hold on
tmp = squeeze(S.Zirep(3, 5, end - nobs + 1 : end));
plot(tmp(1 : lengthToPlot), 'r')


%% daily loadings

dtIrep  = mean(full(dtS), 4);

figure
cnt = 1;
corrZ = nan(1,9);
for i = 1 : 9
   subplot(3, 3, cnt)
   estz = permute(ZquantileAll(3 + i, 7, :, :), [3 4 1 2]);
   plot(estz, 'b')   
   hold on   
   truez = squeeze(S.Zirep(3 + i, 7, end - nobs + 1 : end));
   plot(truez(1 : lengthToPlot), 'r')   
   cnt = cnt + 1;
   
   corrZ(i) = corr(estz(:,2), truez(1 : lengthToPlot));
  
end

mean(corrZ(2:end))

%% Individual loadings to plot

for i = [2 5 8]
        % Individual plots
        figure('name', ['sim_loading_n' num2str(i)]);
        estz = permute(ZquantileAll(3 + i, 7, :, :), [3 4 1 2]);
        plot(estz(:, 2), 'k', 'linewidth', 1)
        hold on
        set(gca, 'ylim', [-2.5, 1.5])
        plot(estz(:, [1 3]), 'k', 'linewidth', 1, 'linestyle','-.')
        plot(dtIrep(3 + 2).*ones(size(estz, 1), 1), 'k')
        plot(-dtIrep(3 + 2).*ones(size(estz, 1), 1), 'k')
        %axis tight
        box off
        axis tight
        set(gca, 'ylim', [-2.5, 1.5])       
        
        figure('name', ['true_loading_n' num2str(i)]);
        truez  = squeeze(S.Zirep(3 + i, 7, end - nobs + 1 : end));                        
        truezu = squeeze(S.Zirepu(3 + i, 7, end - nobs + 1 : end));                        
        p1 = plot(truez(1 : lengthToPlot), 'k', 'linewidth', 1);
        hold on        
        set(gca, 'ylim', [-2.5, 1.5])
        barp = nan.*truez(1 : lengthToPlot);
        barp(truezu(1 : lengthToPlot) == 0) = max(get(gca,'ylim'));
        bsIdx = [];
        for ii = 1 : numel(barp)                        
           if ~isnan(barp(ii))  
               if ii == 1 
                   st = ii;
               else
                   if isnan(barp(ii-1))
                       st = ii;
                   end
               end
               if ii == numel(barp)                        
                   en = ii;
                   bsIdx = cat(1, bsIdx, [st en]);
               else
                   if isnan(barp(ii+1))                  
                       en = ii;
                       bsIdx = cat(1, bsIdx, [st en]);
                   end               
               end
           end
        end        
        shade(bsIdx(:,1), bsIdx(:,2), [0.85 0.85 0.85]);               
        p1 = plot(truez(1 : lengthToPlot), 'k', 'linewidth', 1);        
        plot(0.2.*ones(size(truez(1 : lengthToPlot), 1), 1), 'k')
        plot(-0.2.*ones(size(truez(1 : lengthToPlot), 1), 1), 'k')
        %axis tight
        box off        
        axis tight
        set(gca, 'ylim', [-2.5, 1.5])
        
        
end



%% Individual loadings to plot

Zds     = permute(full(ZS(4:end, 7, :, :)), [3 4 1 2]);
dtIrep  = permute(full(dtS(4:end, 1, 1, :)), [1 4 2 3]);
draws   = size(dtIrep, 2);

for i = [2 5 8]
        % Individual plots
        figure('name', ['sim_thresholdprob_n' num2str(i)]);
        idx = squeeze(...
                    abs(Zds(:, :, i))...
                    < ...
                    repmat(dtIrep(i, :),[nobs 1])...
                    );                             
       Zdsi      = Zds(:, :, i);
       Zdsi(idx) = 0;
       ah = area((sum(Zdsi(1 : lengthToPlot, :) == 0, 2)/draws)', 'facecolor', [0.85 0.85 0.85], 'edgecolor', [0.85 0.85 0.85]);       
       box off        
       axis tight
       set(gca, 'ylim', [0 1]);       
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%