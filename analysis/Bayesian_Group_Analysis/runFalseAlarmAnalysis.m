% 29-07-20 DL
% Figure C1
% 7: S/S to S/L, S/H --> These are all "same" trials, so we can't compute HITS to FALSE ALARMS as 
%       there wouldn't be any HITS. We expect CORRECT REJECTIONS to decrease (i.e., 
%       FALSE ALARMS to increase) with the strength of the bottom. We can just compare the error rate across items
clear all
clc
close all
checkchains = true;

subs = 1:6; 
conditions = [1 1 1 2 2 2];

%% Set Up Data
HR = [];
FA = [];
hits = []; falarms = [];
Ns = []; Nn = [];  
for sidx = 1:numel(subs)
    [subHR, subFA, cols] = extractHandFA(subs(sidx));
    HR = [HR; ones(size(subHR,1),1) * subs(sidx), ones(size(subHR,1),1) * conditions(sidx), subHR];
    FA = [FA; ones(size(subFA,1),1) * subs(sidx), ones(size(subFA,1),1) * conditions(sidx), subFA];
    
    hits =  [hits; 
             subs(sidx), conditions(sidx), 1, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))';
             subs(sidx), conditions(sidx), 2, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 2, strcmp(cols, 'Nresp'))';
             subs(sidx), conditions(sidx), 3, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 3, strcmp(cols, 'Nresp'))'];
    
    Ns = [Ns;
          subs(sidx), conditions(sidx), 1, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))';
          subs(sidx), conditions(sidx), 2, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 2, strcmp(cols, 'Ntot'))';
          subs(sidx), conditions(sidx), 3, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 3, strcmp(cols, 'Ntot'))'];
    
    falarms = [falarms; 
        subs(sidx), conditions(sidx), 2, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))';
        subs(sidx), conditions(sidx), 3, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 2, strcmp(cols, 'Nresp'))';
        subs(sidx), conditions(sidx), 4, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 3, strcmp(cols, 'Nresp'))'];
    
    Nn = [Nn; 
        subs(sidx), conditions(sidx), 2, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))';
        subs(sidx), conditions(sidx), 3, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 2, strcmp(cols, 'Ntot'))';
        subs(sidx), conditions(sidx), 4, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 3, strcmp(cols, 'Ntot'))'];
    
end
cols = ['sub', 'cond', cols]; % 'sub'    'cond'    'top'    'bot'    'dline'    'rate'    'Nresp'    'Ntot'    'rt'

deadlines = [.05, .1, .2, .4, .8, 1.8];
rts = aggregate([HR; FA], strcmp(cols, 'dline'), strcmp(cols, 'rt'));
tpt = deadlines + rts(:,2)'./1000;
D = numel(deadlines); % Number of deadlines
C = numel(unique(hits(:,2))); % Number of conditions
I = numel(unique(hits(:,3))); % Number of item conditions
R = size(hits, 1); % Number of rows (i.e., subjects x conditions x items)

% d[i,j] <- (m * (1 - exp(-a * (t[j] - T0)))) * (1/pow(a * (t[j] - T0) * (pow(s,2) + (1/(t[j] - T0))), .5))

%% MCMC Settings for JAGS
nchains  = 3;     % How Many Chains?
nburnin  = 5000;   % How Many Burn-in Samples?
nsamples = 5000;  % How Many Recorded Samples? originally 5000

%% Assign Matlab Variables to the Observed JAGS Nodes
datastruct = struct('hits', hits(:,4:end), 'falarms', falarms(:,4:end),...
    'Ns', Ns(:,4:end), 'Nn', Nn(:,4:end),...
    'cidx', hits(:,2), 'iidx', hits(:,3),...
    't', tpt, 'R', R, 'D', D, 'C', C, 'I', I);

%% Initialize Unobserved Variables
% Starting points for the unobserved variables in JAGS
for i=1:nchains
    S.mu_m = ones(C, I);
    S.precm = 10 * ones(C,I);
    S.mu_T0 = .05 * ones(C,I);
    S.precT0 = 10 * ones(C,I);
    S.mu_tau = ones(C,I);
    S.prectau = 10 * ones(C,I);
    
    S.m = rand(R,1);
    S.T0 = rand(R,1);
    S.tau = rand(R,1);
    
    S.c = rand(R,6);
    
    init0(i) = S;
end

%% Run JAGS
savefn = 'samples_falseAlarms.mat';
jagsModelFileName = 'stopSignalModel.txt';
whichParmsToMonitor = {'mu_m', 'prec_m', 'mu_T0', 'prec_T0', 'mu_tau', 'prec_tau',...
                       'm', 'T0', 'tau', 'd', 'c', 'dprime',...
                       'msamp', 'T0samp', 'tausamp', 'thetaf'};

if exist(savefn, 'file') ~= 2
    % Pass all information to JAGS
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, jagsModelFileName), ...
        init0, ...
        'doparallel' , false, ...
        'nchains', nchains,...
        'nburnin', nburnin,...
        'nsamples', nsamples, ...
        'thin', 10, ...
        'monitorparams', whichParmsToMonitor, ...
        'savejagsoutput' , 1 , ...
        'verbosity' , 1 , ...
        'cleanup' , 0 , ...
        'workingdir' , 'tmpjags' );
    
    save(savefn, 'samples', 'stats')
else
    load(savefn)
end

%% Check chains for convergence
if checkchains
% fieldnames(stats.mean); 
plotChains(samples.mu_m,...
    makeVarName(samples.mu_m, 'mu m'))
plotChains(samples.prec_m,...
    makeVarName(samples.prec_m, 'prec m'))

plotChains(samples.mu_T0,...
    makeVarName(samples.mu_T0, 'mu T0'))
plotChains(samples.prec_T0,...
    makeVarName(samples.prec_T0, 'prec T0'))

plotChains(samples.mu_tau,...
    makeVarName(samples.mu_tau, 'mu tau'))
plotChains(samples.prec_tau,...
    makeVarName(samples.prec_tau, 'prec tau'))

plotChains(samples.m,...
    makeVarName(samples.m, 'm'))

plotChains(samples.T0,...
    makeVarName(samples.T0, 'T0'))

plotChains(samples.tau,...
    makeVarName(samples.tau, 'tau'))

plotChains(samples.d,...
    makeVarName(samples.d, 'dprime'))
plotChains(samples.c,...
    makeVarName(samples.c, 'c'))
end

% %% Plot posteriors for hyperparms
% % Reshape the samples to combine chains (only works for 3 chains)
% resh = @(x)([samples.(sprintf('%s', x))(1,:,1,1), samples.(sprintf('%s', x))(2,:,1,1), samples.(sprintf('%s', x))(3,:,1,1); % Condition 1, item 1
%              samples.(sprintf('%s', x))(1,:,1,2), samples.(sprintf('%s', x))(2,:,1,2), samples.(sprintf('%s', x))(3,:,1,2); % Condition 1, item 2
%              samples.(sprintf('%s', x))(1,:,1,3), samples.(sprintf('%s', x))(2,:,1,3), samples.(sprintf('%s', x))(3,:,1,3); % Condition 1, item 1
%              samples.(sprintf('%s', x))(1,:,1,4), samples.(sprintf('%s', x))(2,:,1,4), samples.(sprintf('%s', x))(3,:,1,4); % Condition 1, item 2
%              samples.(sprintf('%s', x))(1,:,2,1), samples.(sprintf('%s', x))(2,:,2,1), samples.(sprintf('%s', x))(3,:,2,1); % Condition 1, item 1
%              samples.(sprintf('%s', x))(1,:,2,2), samples.(sprintf('%s', x))(2,:,2,2), samples.(sprintf('%s', x))(3,:,2,2); % Condition 1, item 2
%              samples.(sprintf('%s', x))(1,:,2,3), samples.(sprintf('%s', x))(2,:,2,3), samples.(sprintf('%s', x))(3,:,2,3); % Condition 1, item 1
%              samples.(sprintf('%s', x))(1,:,2,4), samples.(sprintf('%s', x))(2,:,2,4), samples.(sprintf('%s', x))(3,:,2,4); % Condition 1, item 2             
%              ]');  
% m = resh('mu_m'); % Nsamples x (Ncondition x Nitems)
% T0 = resh('mu_T0');
% tau = resh('mu_tau');
% 
% figure('WindowStyle', 'docked')
% colours = [.10 .72 .88; .27 .47 .77; .43 .27 .66; .56 .03 .52;
%            .44 .99 .69; .54 .76 .60; .67 .53 .54; .87 .25 .45];
%        
%          %[.56 .03 .52; .43 .27 .66; .27 .47 .77; .10 .72 .88;
%          % .96 .03 .36; .87 .25 .45; .67 .70 .54; .54 .90 .60];
%            %.56    1 .52; 1 .47 .77;  .27  1 .66; 1 .72 .88];
% 
% subplot(2,6,7); 
% xi = linspace(0, 5, 200);
% for i = 1:4
%     mf(i,:) = ksdensity(m(:,i), xi);
%     hm(i) = fill([xi, fliplr(xi)], [mf(i,:), zeros(1,size(mf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on  
% end
% xlabel('Asymptotic D-Prime')
% ylabel('Posterior Density')
% % legend(hm(1:4), 'Aligned: H/H', 'Aligned: H/L', 'Aligned: L/H', 'Aligned: L/L',...
% %                 'Location', 'Best')
% 
% subplot(2,6,10); 
% for i = 5:size(m,2)
%     mf(i,:) = ksdensity(m(:,i), xi);
%     hm(i) = fill([xi, fliplr(xi)], [mf(i,:), zeros(1,size(mf(i,:),2))], colours(i,:), 'FaceAlpha', .7); hold on  
% end
% xlabel('Asymptotic D-Prime')
% ylabel('Posterior Density')
% % legend(hm(5:end), 'Misaligned: H/H','Misaligned: H/L', 'Misaligned: L/H', 'Misaligned: L/L', 'Location', 'Best')       
%        
% subplot(2,6,8)
% xi = linspace(0, .5, 200);
% for i = 1:4
%     T0f(i,:) = ksdensity(T0(:,i), xi);
%     hT0(i) = fill([xi, fliplr(xi)], [T0f(i,:), zeros(1,size(T0f(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
% end
% xlabel('T0')
% ylabel('Posterior Density')
% %legend(hT0([2 4 1 3]), 'Aligned: Congruent High', 'Aligned: Incongruent High', 'Aligned: Congruent Low', 'Aligned: Incongruent Low')
%        
% subplot(2,6,11)
% for i = 5:8
%     T0f(i,:) = ksdensity(T0(:,i), xi);
%     hT0(i) = fill([xi, fliplr(xi)], [T0f(i,:), zeros(1,size(T0f(i,:),2))], colours(i,:), 'FaceAlpha', .7); hold on
% end
% xlabel('T0')
% ylabel('Posterior Density')
% %legend(hT0([6 8 5 7]), 'Misaligned: Congruent High', 'Misaligned: Incongruent High', 'Misaligned: Congruent Low', 'Misaligned: Incongruent Low')
% 
% subplot(2,6,9)
% xi = linspace(0, 1, 200);
% for i = 1:4
%     tauf(i,:) = ksdensity(tau(:,i), xi);
%     htau(i) = fill([xi, fliplr(xi)], [tauf(i,:), zeros(1,size(tauf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
% end
% xlabel('\tau')
% ylabel('Posterior Density')
% legend(htau([4 3 2 1]), 'C-H', 'I-H', 'C-L', 'I-L')
% % legend(htau, 'Aligned: H/H', 'Aligned: H/L', 'Aligned: L/H', 'Aligned: L/L')
%        
% 
% subplot(2,6,12)
% for i = 5:8
%     tauf(i,:) = ksdensity(tau(:,i), xi);
%     htau(i) = fill([xi, fliplr(xi)], [tauf(i,:), zeros(1,size(tauf(i,:),2))], colours(i,:), 'FaceAlpha', .7); hold on
% end
% xlabel('\tau')
% ylabel('Posterior Density')
% legend(htau([8 7 6 5]), 'C-H', 'I-H', 'C-L', 'I-L')
% % legend(htau, 'Misaligned: H/H','Misaligned: H/L', 'Misaligned: L/H', 'Misaligned: L/L')

%% Use group level parameters to generate group level SAT function predictions
figure('WindowStyle', 'docked');

% Plot SAT function in each condition
fexp = @(m, T0, tau, t)((m .* (1 - exp(-1./tau .* (t - T0)))));
 
xlims = [0 7];
xticks = 1:6;
xticklabs = {'.32','\newline{.36}', '.44', '.64', '1.03', '2.02'};
data_dur = tpt;
t = 0:.001:2.5;


% Map to x-axis space
% xt = [linspace(xlims(1), xticks(1), find(data_dur(1) <= t, 1, 'first')),...                                       % to .05
%     linspace(xticks(1), xticks(2), find(data_dur(2) <= t, 1, 'first')- find(data_dur(1) <= t, 1, 'first')),...  % to .1
%     linspace(xticks(2), xticks(3), find(data_dur(3) <= t, 1, 'first')- find(data_dur(2) <= t, 1, 'first')),...  % to .2
%     linspace(xticks(3), xticks(4), find(data_dur(4) <= t, 1, 'first')- find(data_dur(3) <= t, 1, 'first')),...  % to .4
%     linspace(xticks(4), xticks(5), find(data_dur(5) <= t, 1, 'first')- find(data_dur(4) <= t, 1, 'first')),...  % to .8
%     linspace(xticks(5), xticks(6), find(data_dur(6) <= t, 1, 'first')- find(data_dur(5) <= t, 1, 'first')),...  % to 1.8
%     linspace(xticks(6), xlims(2), find(max(t) <= t, 1, 'first') - find(data_dur(6) <= t, 1, 'first'))];         % to xlim
xt = t; 
xloc = tpt; 
xmax = 2.5;

% nPostSamps = 1000;
% mset = datasample(m, nPostSamps);
% T0set = datasample(T0,nPostSamps);
% tauset = datasample(tau, nPostSamps);
% 
% y11 = cell2mat(arrayfun(@(i)fexp(mset(i,1), T0set(i,1), tauset(i,1), t), 1:nPostSamps, 'uni', false)');
% y12 = cell2mat(arrayfun(@(i)fexp(mset(i,2), T0set(i,2), tauset(i,2), t), 1:nPostSamps, 'uni', false)');
% y13 = cell2mat(arrayfun(@(i)fexp(mset(i,3), T0set(i,3), tauset(i,3), t), 1:nPostSamps, 'uni', false)');
% y14 = cell2mat(arrayfun(@(i)fexp(mset(i,4), T0set(i,4), tauset(i,4), t), 1:nPostSamps, 'uni', false)');
% y21 = cell2mat(arrayfun(@(i)fexp(mset(i,5), T0set(i,5), tauset(i,5), t), 1:nPostSamps, 'uni', false)');
% y22 = cell2mat(arrayfun(@(i)fexp(mset(i,6), T0set(i,6), tauset(i,6), t), 1:nPostSamps, 'uni', false)');
% y23 = cell2mat(arrayfun(@(i)fexp(mset(i,7), T0set(i,7), tauset(i,7), t), 1:nPostSamps, 'uni', false)');
% y24 = cell2mat(arrayfun(@(i)fexp(mset(i,8), T0set(i,8), tauset(i,8), t), 1:nPostSamps, 'uni', false)');

% low11 = prctile(y11, 5); hi11  = prctile(y11, 95);
% low12 = prctile(y12, 5); hi12  = prctile(y12, 95);
% low13 = prctile(y13, 5); hi13  = prctile(y13, 95);
% low14 = prctile(y14, 5); hi14  = prctile(y14, 95);
% low21 = prctile(y21, 5); hi21  = prctile(y21, 95);
% low22 = prctile(y22, 5); hi22  = prctile(y22, 95);
% low23 = prctile(y23, 5); hi23  = prctile(y23, 95);
% low24 = prctile(y24, 5); hi24  = prctile(y24, 95);


% post_d11 = datasample([squeeze(samples.dprime(1,:,:,1,1)); squeeze(samples.dprime(2,:,:,1,1)); squeeze(samples.dprime(3,:,:,1,1))], 500);
% post_d12 = datasample([squeeze(samples.dprime(1,:,:,1,2)); squeeze(samples.dprime(2,:,:,1,2)); squeeze(samples.dprime(3,:,:,1,2))], 500);
% post_d21 = datasample([squeeze(samples.dprime(1,:,:,2,1)); squeeze(samples.dprime(2,:,:,2,1)); squeeze(samples.dprime(3,:,:,2,1))], 500);
% post_d22 = datasample([squeeze(samples.dprime(1,:,:,2,2)); squeeze(samples.dprime(2,:,:,2,2)); squeeze(samples.dprime(3,:,:,2,2))], 500);
post_d11 = datasample([squeeze(samples.thetaf(1,:,1,:)); squeeze(samples.thetaf(1,:,4,:)); squeeze(samples.thetaf(1,:,7,:))], 500);
post_d12 = datasample([squeeze(samples.thetaf(1,:,2,:)); squeeze(samples.thetaf(1,:,5,:)); squeeze(samples.thetaf(1,:,8,:))], 500);
post_d13 = datasample([squeeze(samples.thetaf(1,:,3,:)); squeeze(samples.thetaf(1,:,6,:)); squeeze(samples.thetaf(1,:,9,:))], 500);

post_d21 = datasample([squeeze(samples.thetaf(1,:,10,:)); squeeze(samples.thetaf(1,:,13,:)); squeeze(samples.thetaf(1,:,16,:))], 500);
post_d22 = datasample([squeeze(samples.thetaf(1,:,11,:)); squeeze(samples.thetaf(1,:,14,:)); squeeze(samples.thetaf(1,:,17,:))], 500);
post_d23 = datasample([squeeze(samples.thetaf(1,:,12,:)); squeeze(samples.thetaf(1,:,15,:)); squeeze(samples.thetaf(1,:,18,:))], 500);
post_d11(post_d11==0) = nan;
post_d12(post_d12==0) = nan;
post_d13(post_d13==0) = nan;
post_d21(post_d21==0) = nan;
post_d22(post_d22==0) = nan;
post_d23(post_d23==0) = nan;

subplot(2,2,1)
vh11 = violin(post_d11, 'x', xloc, ...
    'facecolor', [.27 .47 .77], 'edgecolor', 'k', 'facealpha', .7); hold on %C-L
vh12 = violin(post_d12, 'x', xloc, ...
    'facecolor', [.43 .27 .66], 'edgecolor', 'k', 'facealpha', .7); %C-H
vh13 = violin(post_d13, 'x', xloc, ...
    'facecolor', [.56 .03 .52], 'edgecolor', 'k', 'facealpha', .7); %I-L

set(gca, 'XLim', [0 2.5], 'YLim', [0 0.6], 'XTick', 0:.5:2.5, 'XTickLabel',  {'.00', '.50', '1.0', '1.5', '2.0', '2.5'})
% set(gca, 'XLim', [0 xmax], 'YLim', [0 5], 'TickLabelInterpreter', 'tex', 'XTick', xloc, 'XTickLabel',  xticklabs) % num2str(round(tpt,2)')
% set(gca, 'Children', [vh11, vh12, fh11, fh12])
xlabel('Total Processing Time (sec)')
ylabel('False Alarm Rate')
legend([vh13(1), vh12(1), vh11(1)],...
    'Aligned: Same/High','Aligned: Same/Low', 'Aligned: Same/Same',...
    'Location', 'Best')

subplot(2,2,2)
% fh21 = fill([xt, fliplr(xt)], [low21, fliplr(hi21)], [1 0 0], 'FaceAlpha', .2, 'EdgeAlpha', 0); hold on
% fh22 = fill([xt, fliplr(xt)], [low22, fliplr(hi22)], [0 0 1], 'FaceAlpha', .2, 'EdgeAlpha', 0);
% fh23 = fill([xt, fliplr(xt)], [low23, fliplr(hi23)], [1 0 1], 'FaceAlpha', .2, 'EdgeAlpha', 0);
% fh24 = fill([xt, fliplr(xt)], [low24, fliplr(hi24)], [0 1 1], 'FaceAlpha', .2, 'EdgeAlpha', 0);
% .44 .99 .69; .54 .76 .60; .67 .53 .54; .87 .25 .45

vh21 = violin(post_d21, 'x', xloc, ...
    'facecolor', [.54 .76 .60], 'edgecolor', 'k', 'facealpha', .7); hold on
vh22 = violin(post_d22, 'x', xloc, ...
    'facecolor', [.67 .53 .54], 'edgecolor', 'k', 'facealpha', .7);
vh23 = violin(post_d23, 'x', xloc, ...
    'facecolor', [.87 .25 .45], 'edgecolor', 'k', 'facealpha', .7);
% vh24 = violin(post_d24, 'x', xloc, ...
%     'facecolor', [1 1 1], 'edgecolor', [0 .5 .5], 'facealpha', 1);

set(gca, 'XLim', [0 2.5], 'YLim', [0 0.6], 'XTick', 0:.5:2.5, 'XTickLabel',  {'.00', '.50', '1.0', '1.5', '2.0', '2.5'})
% set(gca, 'XLim', [0 xmax], 'YLim', [0 5],'TickLabelInterpreter', 'tex', 'XTick', xloc, 'XTickLabel',  xticklabs)
% set(gca, 'Children', [vh21, vh22, fh21, fh22])
xlabel('Total Processing Time (sec)')
ylabel('False Alarm Rate')
legend([vh23(1), vh22(1), vh21(1)],...
    'Misaligned: Same/High','Misaligned: Same/Low', 'Misaligned: Same/Same',...
    'Location', 'Best')

%% Plot individual table
subjectIds = {'1', '2', '3', '1', '2', '3'};
alignmentConditionIds = {'   Aligned', 'Misaligned'};
itemConditionIds = {'', 'Same/Same', ' Same/Low', 'Same/High'};
totalProcessingTimeStr = cellstr(num2str(round(tpt * 1000, 0)'))';

cons  = falarms(:,1:3);
means = stats.mean.thetaf;
ci_lo = reshape(stats.ci_low.thetaf, R, D);
ci_hi = reshape(stats.ci_high.thetaf, R, D);

cnt = 1; 
for i = 1:R
    table{i,1} = sprintf('%s\t %s\t %4s\t', subjectIds{cons(i,1)}, alignmentConditionIds{cons(i,2)}, itemConditionIds{cons(i,3)});
    for j = 1:D
        table{i,1} = [table{i,1}, sprintf(' .%02d [.%02d-.%02d]\t', round(100*means(i,j),0), round(100*ci_lo(i,j),0), round(100*ci_hi(i,j), 0))];
    end
end