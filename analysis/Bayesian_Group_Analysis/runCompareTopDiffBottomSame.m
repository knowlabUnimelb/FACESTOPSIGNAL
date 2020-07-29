% 29-07-20 DL
% Figure D1
% 1:  S/S to L/S, H/S --> What is the effect on sensitivity of changing the 
%       top face while holding the bottom constant and unchanged?

clear all
clc
close all

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
             subs(sidx), conditions(sidx), 1, subHR(subHR(:,strcmp(cols, 'top')) == 2 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))';
             subs(sidx), conditions(sidx), 2, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))'];
    
    Ns = [Ns;
          subs(sidx), conditions(sidx), 1, subHR(subHR(:,strcmp(cols, 'top')) == 2 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))';
          subs(sidx), conditions(sidx), 2, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))'];
    
    falarms = [falarms; 
        subs(sidx), conditions(sidx), 1, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))';
        subs(sidx), conditions(sidx), 2, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))'];
    
    Nn = [Nn; 
        subs(sidx), conditions(sidx), 1, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))';
        subs(sidx), conditions(sidx), 2, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))'];
    
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
savefn = 'samples_topDiffBottomSame.mat';
jagsModelFileName = 'stopSignalModel.txt';
whichParmsToMonitor = {'mu_m', 'prec_m', 'mu_T0', 'prec_T0', 'mu_tau', 'prec_tau',...
                       'm', 'T0', 'tau', 'd', 'c', 'dprime',...
                       'msamp', 'T0samp', 'tausamp'};

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

%% Plot posteriors for hyperparms
% Reshape the samples to combine chains (only works for 3 chains)
resh = @(x)([samples.(sprintf('%s', x))(1,:,1,1), samples.(sprintf('%s', x))(2,:,1,1), samples.(sprintf('%s', x))(3,:,1,1); % Condition 1, item 1
             samples.(sprintf('%s', x))(1,:,1,2), samples.(sprintf('%s', x))(2,:,1,2), samples.(sprintf('%s', x))(3,:,1,2); % Condition 1, item 2
             samples.(sprintf('%s', x))(1,:,2,1), samples.(sprintf('%s', x))(2,:,2,1), samples.(sprintf('%s', x))(3,:,2,1); % Condition 1, item 1
             samples.(sprintf('%s', x))(1,:,2,2), samples.(sprintf('%s', x))(2,:,2,2), samples.(sprintf('%s', x))(3,:,2,2); % Condition 1, item 2
             ]');  
m = resh('mu_m'); % Nsamples x (Ncondition x Nitems)
T0 = resh('mu_T0');
tau = resh('mu_tau');

figure('WindowStyle', 'docked')
colours = [.27 .47 .77; .56 .03 .52; .54 .76 .60; .87 .25 .45];
         % [.10 .72 .88; .43 .27 .66; .44 .99 .69; .67 .53 .54];   
         
%subplot(2,3,4)
%xi = linspace(0, 5, 200);
%for i = 1:size(m,2)
 %   mf(i,:) = ksdensity(m(:,i), xi);
  %  hm(i) = fill([xi, fliplr(xi)], [mf(i,:), zeros(1,size(mf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on  
%end
%xlabel('Asymptotic D-Prime')
%ylabel('Posterior Density')
%legend(hm, 'Aligned: L/S', 'Aligned: H/S','Misaligned: L/S', 'Misaligned: H/S')

%subplot(2,3,5)
%xi = linspace(0, .5, 200);
%for i = 1:size(T0,2)
%    T0f(i,:) = ksdensity(T0(:,i), xi);
 %   hT0(i) = fill([xi, fliplr(xi)], [T0f(i,:), zeros(1,size(T0f(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
%end
%xlabel('T0')
%ylabel('Posterior Density')
%legend(hT0, 'Aligned: L/S', 'Aligned: H/S','Misaligned: L/S', 'Misaligned: H/S')

%subplot(2,3,6)
%xi = linspace(0, 1, 200);
%for i = 1:size(tau,2)
 %   tauf(i,:) = ksdensity(tau(:,i), xi);
  %  htau(i) = fill([xi, fliplr(xi)], [tauf(i,:), zeros(1,size(tauf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
%end
%xlabel('\tau')
%ylabel('Posterior Density')
%legend(htau, 'Aligned: L/S', 'Aligned: H/S','Misaligned: L/S', 'Misaligned: H/S')

subplot(2,6,7); 
xi = linspace(0, 5, 200);
for i = 1:2 %1:4
    mf(i,:) = ksdensity(m(:,i), xi);
    hm(i) = fill([xi, fliplr(xi)], [mf(i,:), zeros(1,size(mf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on  
end
xlabel('Asymptotic D-Prime')
ylabel('Posterior Density')
% legend(hm(1:4), 'Aligned: H/H', 'Aligned: H/L', 'Aligned: L/H', 'Aligned: L/L',...
%                 'Location', 'Best')

subplot(2,6,10); 
for i = 3:size(m,2) %5:size(m,2)
    mf(i,:) = ksdensity(m(:,i), xi);
    hm(i) = fill([xi, fliplr(xi)], [mf(i,:), zeros(1,size(mf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on  
end
xlabel('Asymptotic D-Prime')
ylabel('Posterior Density')
% legend(hm(5:end), 'Misaligned: H/H','Misaligned: H/L', 'Misaligned: L/H', 'Misaligned: L/L', 'Location', 'Best')       
       
subplot(2,6,8)
xi = linspace(0, .5, 200);
for i = 1:2 %1:4
    T0f(i,:) = ksdensity(T0(:,i), xi);
    hT0(i) = fill([xi, fliplr(xi)], [T0f(i,:), zeros(1,size(T0f(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
end
xlabel('T0')
ylabel('Posterior Density')
%legend(hT0(1:4), 'Aligned: L/L', 'Aligned: L/H', 'Aligned: H/L', 'Aligned: H/H')
       
subplot(2,6,11)
for i = 3:4 %5:8
    T0f(i,:) = ksdensity(T0(:,i), xi);
    hT0(i) = fill([xi, fliplr(xi)], [T0f(i,:), zeros(1,size(T0f(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
end
xlabel('T0')
ylabel('Posterior Density')
%legend(hT0(5:end), 'Misaligned: L/L','Misaligned: L/H', 'Misaligned: H/L', 'Misaligned: H/H')
       

subplot(2,6,9)
xi = linspace(0, 1, 200);
for i = 1:2 %1:4
    tauf(i,:) = ksdensity(tau(:,i), xi);
    htau(i) = fill([xi, fliplr(xi)], [tauf(i,:), zeros(1,size(tauf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
end
xlabel('\tau')
ylabel('Posterior Density')
legend(htau([2 1]), 'H/S', 'L/S')
%legend(htau([4 3 2 1]), 'H/H', 'H/L', 'L/H', 'L/L')
       

subplot(2,6,12)
for i = 3:4 %5:8
    tauf(i,:) = ksdensity(tau(:,i), xi);
    htau(i) = fill([xi, fliplr(xi)], [tauf(i,:), zeros(1,size(tauf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
end
xlabel('\tau')
ylabel('Posterior Density')
legend(htau([4 3]), 'H/S','L/S')
% legend(htau([8 7 6 5]), 'H/H','H/L', 'L/H', 'L/L')

%% Use group level parameters to generate group level SAT function predictions
% figure('WindowStyle', 'docked');

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

nPostSamps = 1000;
mset = datasample(m, nPostSamps);
T0set = datasample(T0,nPostSamps);
tauset = datasample(tau, nPostSamps);

y11 = cell2mat(arrayfun(@(i)fexp(mset(i,1), T0set(i,1), tauset(i,1), t), 1:nPostSamps, 'uni', false)');
y12 = cell2mat(arrayfun(@(i)fexp(mset(i,2), T0set(i,2), tauset(i,2), t), 1:nPostSamps, 'uni', false)');
y21 = cell2mat(arrayfun(@(i)fexp(mset(i,3), T0set(i,3), tauset(i,3), t), 1:nPostSamps, 'uni', false)');
y22 = cell2mat(arrayfun(@(i)fexp(mset(i,4), T0set(i,4), tauset(i,4), t), 1:nPostSamps, 'uni', false)');

low11 = prctile(y11, 5); hi11  = prctile(y11, 95);
low12 = prctile(y12, 5); hi12  = prctile(y12, 95);
low21 = prctile(y21, 5); hi21  = prctile(y21, 95);
low22 = prctile(y22, 5); hi22  = prctile(y22, 95);



% post_d11 = datasample([squeeze(samples.dprime(1,:,:,1,1)); squeeze(samples.dprime(2,:,:,1,1)); squeeze(samples.dprime(3,:,:,1,1))], 500);
% post_d12 = datasample([squeeze(samples.dprime(1,:,:,1,2)); squeeze(samples.dprime(2,:,:,1,2)); squeeze(samples.dprime(3,:,:,1,2))], 500);
% post_d21 = datasample([squeeze(samples.dprime(1,:,:,2,1)); squeeze(samples.dprime(2,:,:,2,1)); squeeze(samples.dprime(3,:,:,2,1))], 500);
% post_d22 = datasample([squeeze(samples.dprime(1,:,:,2,2)); squeeze(samples.dprime(2,:,:,2,2)); squeeze(samples.dprime(3,:,:,2,2))], 500);
post_d11 = datasample([squeeze(samples.d(1,:,1,:)); squeeze(samples.d(1,:,3,:)); squeeze(samples.d(1,:,5,:))], 500);
post_d12 = datasample([squeeze(samples.d(1,:,2,:)); squeeze(samples.d(1,:,4,:)); squeeze(samples.d(1,:,6,:))], 500);
post_d21 = datasample([squeeze(samples.d(1,:,7,:)); squeeze(samples.d(1,:,9,:)); squeeze(samples.d(1,:,11,:))], 500);
post_d22 = datasample([squeeze(samples.d(1,:,8,:)); squeeze(samples.d(1,:,10,:)); squeeze(samples.d(1,:,12,:))], 500);

post_d11(post_d11==0) = nan;
post_d12(post_d12==0) = nan;
post_d21(post_d21==0) = nan;
post_d22(post_d22==0) = nan;

subplot(2,2,1)
fh11 = fill([xt, fliplr(xt)], [low11, fliplr(hi11)], [.27 .47 .77], 'FaceAlpha', .5, 'EdgeAlpha', 0); hold on
fh12 = fill([xt, fliplr(xt)], [low12, fliplr(hi12)], [.56 .03 .52], 'FaceAlpha', .5, 'EdgeAlpha', 0);
vh11 = violin(post_d11, 'x', xloc, ...
    'facecolor', [.27 .47 .77], 'edgecolor', 'k', 'facealpha', .8); hold on
vh12 = violin(post_d12, 'x', xloc, ...
    'facecolor', [.56 .03 .52], 'edgecolor', 'k', 'facealpha', .8);
set(gca, 'XLim', [0 2.5], 'YLim', [0 5], 'XTick', 0:.5:2.5, 'XTickLabel',  {'.00', '.50', '1.0', '1.5', '2.0', '2.5'})
% set(gca, 'XLim', [0 xmax], 'YLim', [0 5], 'TickLabelInterpreter', 'tex', 'XTick', xloc, 'XTickLabel',  xticklabs) % num2str(round(tpt,2)')
% set(gca, 'Children', [vh11, vh12, fh11, fh12])
xlabel('Total Processing Time (sec)')
ylabel('d-prime')
legend([fh12, fh11], 'Aligned: High/Same', 'Aligned: Low/Same', 'Location', 'Best')

subplot(2,2,2)
fh21 = fill([xt, fliplr(xt)], [low21, fliplr(hi21)], [.54 .76 .60], 'FaceAlpha', .5, 'EdgeAlpha', 0); hold on
fh22 = fill([xt, fliplr(xt)], [low22, fliplr(hi22)], [.87 .25 .45], 'FaceAlpha', .5, 'EdgeAlpha', 0);
vh21 = violin(post_d21, 'x', xloc, ...
    'facecolor', [.54 .76 .60], 'edgecolor', 'k', 'facealpha', .8); hold on
vh22 = violin(post_d22, 'x', xloc, ...
    'facecolor', [.87 .25 .45], 'edgecolor', 'k', 'facealpha', .8);
set(gca, 'XLim', [0 2.5], 'YLim', [0 5], 'XTick', 0:.5:2.5, 'XTickLabel',  {'.00', '.50', '1.0', '1.5', '2.0', '2.5'})
% set(gca, 'XLim', [0 xmax], 'YLim', [0 5],'TickLabelInterpreter', 'tex', 'XTick', xloc, 'XTickLabel',  xticklabs)
% set(gca, 'Children', [vh21, vh22, fh21, fh22])
xlabel('Total Processing Time (sec)')
ylabel('d-prime')
legend([fh22, fh21], 'Misaligned: High/Same', 'Misaligned: Low/Same', 'Location', 'Best')