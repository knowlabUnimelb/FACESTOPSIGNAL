% 29-07-20 DL
% comp8 - Figures 5, 6, and 7
% congruency x alignment interaction

clear all
clc
close all
checkchains = true;

subs = 1:6; 
conditions = [1 1 1 2 2 2];
% 
% subs = [1 4];
% subs = [2 5];
% subs = [3 6];
% conditions = [1 2];
% 
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
             subs(sidx), conditions(sidx), 1, subHR(subHR(:,strcmp(cols, 'top')) == 2 & subHR(:,strcmp(cols, 'bot')) == 2, strcmp(cols, 'Nresp'))';
             subs(sidx), conditions(sidx), 2, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 3, strcmp(cols, 'Nresp'))';
             subs(sidx), conditions(sidx), 3, subHR(subHR(:,strcmp(cols, 'top')) == 2 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))';
             subs(sidx), conditions(sidx), 4, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))'];
    
    Ns = [Ns;
          subs(sidx), conditions(sidx), 1, subHR(subHR(:,strcmp(cols, 'top')) == 2 & subHR(:,strcmp(cols, 'bot')) == 2, strcmp(cols, 'Ntot'))';
          subs(sidx), conditions(sidx), 2, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 3, strcmp(cols, 'Ntot'))';
          subs(sidx), conditions(sidx), 3, subHR(subHR(:,strcmp(cols, 'top')) == 2 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))';
          subs(sidx), conditions(sidx), 4, subHR(subHR(:,strcmp(cols, 'top')) == 3 & subHR(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))'];
    
    falarms = [falarms; 
        subs(sidx), conditions(sidx), 1, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))';
        subs(sidx), conditions(sidx), 2, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Nresp'))';
        subs(sidx), conditions(sidx), 3, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 2, strcmp(cols, 'Nresp'))';
        subs(sidx), conditions(sidx), 4, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 3, strcmp(cols, 'Nresp'))'];
    
    Nn = [Nn; 
        subs(sidx), conditions(sidx), 1, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))';
        subs(sidx), conditions(sidx), 2, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 1, strcmp(cols, 'Ntot'))';
        subs(sidx), conditions(sidx), 3, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 2, strcmp(cols, 'Ntot'))';
        subs(sidx), conditions(sidx), 4, subFA(subFA(:,strcmp(cols, 'top')) == 1 & subFA(:,strcmp(cols, 'bot')) == 3, strcmp(cols, 'Ntot'))'];
    
end
cols = ['sub', 'cond', cols]; % 'sub'    'cond'    'top'    'bot'    'dline'    'rate'    'Nresp'    'Ntot'    'rt'

deadlines = [.05, .1, .2, .4, .8, 1.8];
rts = aggregate([HR; FA], strcmp(cols, 'dline'), strcmp(cols, 'rt'));
tpt = deadlines + rts(:,2)'./1000;
D = numel(deadlines); % Number of deadlines
C = numel(unique(hits(:,2))); % Number of alignment conditions
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
savefn = 'samples_congruencyAlignmentInteraction.mat';
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

%% Plot posteriors for hyperparms
% Reshape the samples to combine chains (only works for 3 chains)
resh = @(x)([samples.(sprintf('%s', x))(1,:,1,1), samples.(sprintf('%s', x))(2,:,1,1), samples.(sprintf('%s', x))(3,:,1,1); % Condition 1, item 1
             samples.(sprintf('%s', x))(1,:,1,2), samples.(sprintf('%s', x))(2,:,1,2), samples.(sprintf('%s', x))(3,:,1,2); % Condition 1, item 2
             samples.(sprintf('%s', x))(1,:,1,3), samples.(sprintf('%s', x))(2,:,1,3), samples.(sprintf('%s', x))(3,:,1,3); % Condition 1, item 1
             samples.(sprintf('%s', x))(1,:,1,4), samples.(sprintf('%s', x))(2,:,1,4), samples.(sprintf('%s', x))(3,:,1,4); % Condition 1, item 2
             samples.(sprintf('%s', x))(1,:,2,1), samples.(sprintf('%s', x))(2,:,2,1), samples.(sprintf('%s', x))(3,:,2,1); % Condition 1, item 1
             samples.(sprintf('%s', x))(1,:,2,2), samples.(sprintf('%s', x))(2,:,2,2), samples.(sprintf('%s', x))(3,:,2,2); % Condition 1, item 2
             samples.(sprintf('%s', x))(1,:,2,3), samples.(sprintf('%s', x))(2,:,2,3), samples.(sprintf('%s', x))(3,:,2,3); % Condition 1, item 1
             samples.(sprintf('%s', x))(1,:,2,4), samples.(sprintf('%s', x))(2,:,2,4), samples.(sprintf('%s', x))(3,:,2,4); % Condition 1, item 2             
             ]');  
m = resh('mu_m'); % Nsamples x (Ncondition x Nitems)
T0 = resh('mu_T0');
tau = resh('mu_tau');

figure('WindowStyle', 'docked')
colours = [.10 .72 .88; .27 .47 .77; .43 .27 .66; .56 .03 .52; 
           .44 .99 .69; .54 .76 .60; .67 .53 .54; .87 .25 .45];

          %[.97 .93 .22; .96 .73 .27; .96 .73 .27; .96 .56 .33;
          % 1 1 .07; .97 .93 .22; .96 .73 .27; .96 .56 .33]

          %[.97 .93 .22; .96 .73 .27; .96 .56 .33; .96 .41 .40;
          % .97 .93 .22; .96 .73 .27; .96 .56 .33; .96 .41 .40]

          % [.10 .72 .88; .27 .47 .77; .43 .27 .66; .56 .03 .52; 
          %  .44 .99 .69; .54 .76 .60; .67 .53 .54; .87 .25 .45]; original colour scheme
       
subplot(2,6,7); 
xi = linspace(0, 5, 200);
for i = 1:4
    mf(i,:) = ksdensity(m(:,i), xi);
    hm(i) = fill([xi, fliplr(xi)], [mf(i,:), zeros(1,size(mf(i,:),2))], colours(i,:), 'FaceAlpha', .7); hold on  
end
xlabel('Asymptotic D-Prime')
ylabel('Posterior Density')
% legend(hm(1:4), 'Aligned: H/H', 'Aligned: H/L', 'Aligned: L/H', 'Aligned: L/L',...
%                 'Location', 'Best')

subplot(2,6,10); 
for i = 5:size(m,2)
    mf(i,:) = ksdensity(m(:,i), xi);
    hm(i) = fill([xi, fliplr(xi)], [mf(i,:), zeros(1,size(mf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on  
end
xlabel('Asymptotic D-Prime')
ylabel('Posterior Density')
% legend(hm(5:end), 'Misaligned: H/H','Misaligned: H/L', 'Misaligned: L/H', 'Misaligned: L/L', 'Location', 'Best')       
       
subplot(2,6,8)
xi = linspace(0, .5, 200);
for i = 1:4
    T0f(i,:) = ksdensity(T0(:,i), xi);
    hT0(i) = fill([xi, fliplr(xi)], [T0f(i,:), zeros(1,size(T0f(i,:),2))], colours(i,:), 'FaceAlpha', .7); hold on
end
xlabel('T0')
ylabel('Posterior Density')
%legend(hT0([2 4 1 3]), 'Aligned: Congruent High', 'Aligned: Incongruent High', 'Aligned: Congruent Low', 'Aligned: Incongruent Low')
       
subplot(2,6,11)
for i = 5:8
    T0f(i,:) = ksdensity(T0(:,i), xi);
    hT0(i) = fill([xi, fliplr(xi)], [T0f(i,:), zeros(1,size(T0f(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
end
xlabel('T0')
ylabel('Posterior Density')
%legend(hT0([6 8 5 7]), 'Misaligned: Congruent High', 'Misaligned: Incongruent High', 'Misaligned: Congruent Low', 'Misaligned: Incongruent Low')

subplot(2,6,9)
xi = linspace(0, 1, 200);
for i = 1:4
    tauf(i,:) = ksdensity(tau(:,i), xi);
    htau(i) = fill([xi, fliplr(xi)], [tauf(i,:), zeros(1,size(tauf(i,:),2))], colours(i,:), 'FaceAlpha', .7); hold on
end
xlabel('\tau')
ylabel('Posterior Density')
legend(htau([1 2 3 4]), 'C-L', 'C-H', 'I-L', 'I-H')
% legend(htau, 'Aligned: H/H', 'Aligned: H/L', 'Aligned: L/H', 'Aligned: L/L')
       

subplot(2,6,12)
for i = 5:8
    tauf(i,:) = ksdensity(tau(:,i), xi);
    htau(i) = fill([xi, fliplr(xi)], [tauf(i,:), zeros(1,size(tauf(i,:),2))], colours(i,:), 'FaceAlpha', .5); hold on
end
xlabel('\tau')
ylabel('Posterior Density')
legend(htau([5 6 7 8]), 'C-L', 'C-H', 'I-L', 'I-H')
% legend(htau, 'Misaligned: H/H','Misaligned: H/L', 'Misaligned: L/H', 'Misaligned: L/L')

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

nPostSamps = 5000; %200000;
mset = datasample(m, nPostSamps);
T0set = datasample(T0,nPostSamps);
tauset = datasample(tau, nPostSamps);

y11 = cell2mat(arrayfun(@(i)fexp(mset(i,1), T0set(i,1), tauset(i,1), t), 1:nPostSamps, 'uni', false)');
y12 = cell2mat(arrayfun(@(i)fexp(mset(i,2), T0set(i,2), tauset(i,2), t), 1:nPostSamps, 'uni', false)');
y13 = cell2mat(arrayfun(@(i)fexp(mset(i,3), T0set(i,3), tauset(i,3), t), 1:nPostSamps, 'uni', false)');
y14 = cell2mat(arrayfun(@(i)fexp(mset(i,4), T0set(i,4), tauset(i,4), t), 1:nPostSamps, 'uni', false)');
y21 = cell2mat(arrayfun(@(i)fexp(mset(i,5), T0set(i,5), tauset(i,5), t), 1:nPostSamps, 'uni', false)');
y22 = cell2mat(arrayfun(@(i)fexp(mset(i,6), T0set(i,6), tauset(i,6), t), 1:nPostSamps, 'uni', false)');
y23 = cell2mat(arrayfun(@(i)fexp(mset(i,7), T0set(i,7), tauset(i,7), t), 1:nPostSamps, 'uni', false)');
y24 = cell2mat(arrayfun(@(i)fexp(mset(i,8), T0set(i,8), tauset(i,8), t), 1:nPostSamps, 'uni', false)');

low11 = prctile(y11, 5); hi11  = prctile(y11, 95);
low12 = prctile(y12, 5); hi12  = prctile(y12, 95);
low13 = prctile(y13, 5); hi13  = prctile(y13, 95);
low14 = prctile(y14, 5); hi14  = prctile(y14, 95);
low21 = prctile(y21, 5); hi21  = prctile(y21, 95);
low22 = prctile(y22, 5); hi22  = prctile(y22, 95);
low23 = prctile(y23, 5); hi23  = prctile(y23, 95);
low24 = prctile(y24, 5); hi24  = prctile(y24, 95);

% To compute the interaction posterior, I need to take a linear combination
% of the d-prime samples for each congruency x alignment interaction
%   sample.d is in order of alignment x subject x item (congruent_low,
%   congruent_high, incongruent_low, incongruent_high) with condition
%   changing the fastest
%
% Dimensions of samples.d are chain, sample, condition (see first 3 columns
% of falarms), and response signal
%
% To find the interaction posterior, I combine the chains and samples from
% each subject and then downsample those. I can then compute the
% interaction and have a set of interaction samples. 

    % Aligned_congruent_low
%     post_d11 = datasample([[squeeze(samples.d(1,:,1,:)); squeeze(samples.d(1,:,5,:)); squeeze(samples.d(1,:,9,:))];
%                            [squeeze(samples.d(2,:,1,:)); squeeze(samples.d(2,:,5,:)); squeeze(samples.d(2,:,9,:))];
%                            [squeeze(samples.d(3,:,1,:)); squeeze(samples.d(3,:,5,:)); squeeze(samples.d(3,:,9,:))]], 1500);

nt = nchains*nsamples;
acl      = [1*ones(nt,1), squeeze([samples.d(1,:,1,:), samples.d(2,:,1,:), samples.d(3,:,1,:)]); % P1
            2*ones(nt,1), squeeze([samples.d(1,:,5,:), samples.d(2,:,5,:), samples.d(3,:,5,:)]); % P2
            3*ones(nt,1), squeeze([samples.d(1,:,9,:), samples.d(2,:,9,:), samples.d(3,:,9,:)])]; % P3
post_dacl = acl(:,2:end); %datasample(acl(:,2:end), 1500);

    % Aligned_congruent_high
%     post_d12 = datasample([[squeeze(samples.d(1,:,2,:)); squeeze(samples.d(1,:,6,:)); squeeze(samples.d(1,:,10,:))];
%         [squeeze(samples.d(2,:,2,:)); squeeze(samples.d(2,:,6,:)); squeeze(samples.d(2,:,10,:))];
%         [squeeze(samples.d(3,:,2,:)); squeeze(samples.d(3,:,6,:)); squeeze(samples.d(3,:,10,:))]], 1500);

ach      = [1*ones(nt,1), squeeze([samples.d(1,:,2,:), samples.d(2,:,2,:), samples.d(3,:,2,:)]); % P1
            2*ones(nt,1), squeeze([samples.d(1,:,6,:), samples.d(2,:,6,:), samples.d(3,:,6,:)]); % P2
            3*ones(nt,1), squeeze([samples.d(1,:,10,:), samples.d(2,:,10,:), samples.d(3,:,10,:)])]; % P3
post_dach = ach(:,2:end); %datasample(ach(:,2:end), 1500);

    % Aligned_incongruent_low
%     post_d13 = datasample([[squeeze(samples.d(1,:,3,:)); squeeze(samples.d(1,:,7,:)); squeeze(samples.d(1,:,11,:))];
%         [squeeze(samples.d(2,:,3,:)); squeeze(samples.d(2,:,7,:)); squeeze(samples.d(2,:,11,:))];
%         [squeeze(samples.d(3,:,3,:)); squeeze(samples.d(3,:,7,:)); squeeze(samples.d(3,:,11,:))]], 1500);

ail      = [1*ones(nt,1), squeeze([samples.d(1,:,3,:), samples.d(2,:,3,:), samples.d(3,:,3,:)]); % P1
            2*ones(nt,1), squeeze([samples.d(1,:,7,:), samples.d(2,:,7,:), samples.d(3,:,7,:)]); % P2
            3*ones(nt,1), squeeze([samples.d(1,:,11,:), samples.d(2,:,11,:), samples.d(3,:,11,:)])]; % P3
post_dail = ail(:,2:end); %datasample(ail(:,2:end), 1500);

    % Aligned_incongruent_high
%     post_d14 = datasample([[squeeze(samples.d(1,:,4,:)); squeeze(samples.d(1,:,8,:)); squeeze(samples.d(1,:,12,:))];
%         [squeeze(samples.d(2,:,4,:)); squeeze(samples.d(2,:,8,:)); squeeze(samples.d(2,:,12,:))];
%         [squeeze(samples.d(3,:,4,:)); squeeze(samples.d(3,:,8,:)); squeeze(samples.d(3,:,12,:))]], 1500);

aih      = [1*ones(nt,1), squeeze([samples.d(1,:,4,:), samples.d(2,:,4,:), samples.d(3,:,4,:)]); % P1
            2*ones(nt,1), squeeze([samples.d(1,:,8,:), samples.d(2,:,8,:), samples.d(3,:,8,:)]); % P2
            3*ones(nt,1), squeeze([samples.d(1,:,12,:), samples.d(2,:,12,:), samples.d(3,:,12,:)])]; % P3
post_daih = aih(:,2:end); %datasample(aih(:,2:end), 1500);
    
    
    % Misaligned_congruent_low    
%     post_d21 = datasample([[squeeze(samples.d(1,:,13,:)); squeeze(samples.d(1,:,17,:)); squeeze(samples.d(1,:,21,:))];
%         [squeeze(samples.d(2,:,13,:)); squeeze(samples.d(2,:,17,:)); squeeze(samples.d(2,:,21,:))];
%         [squeeze(samples.d(3,:,13,:)); squeeze(samples.d(3,:,17,:)); squeeze(samples.d(3,:,21,:))]], 1500);

mcl      = [1*ones(nt,1), squeeze([samples.d(1,:,13,:), samples.d(2,:,13,:), samples.d(3,:,13,:)]); % P1
            2*ones(nt,1), squeeze([samples.d(1,:,17,:), samples.d(2,:,17,:), samples.d(3,:,17,:)]); % P2
            3*ones(nt,1), squeeze([samples.d(1,:,21,:), samples.d(2,:,21,:), samples.d(3,:,21,:)])]; % P3
post_dmcl = mcl(:,2:end); %datasample(mcl(:,2:end), 1500);


    % Misaligned_congruent_high        
%     post_d22 = datasample([[squeeze(samples.d(1,:,14,:)); squeeze(samples.d(1,:,18,:)); squeeze(samples.d(1,:,22,:))];
%         [squeeze(samples.d(2,:,14,:)); squeeze(samples.d(2,:,18,:)); squeeze(samples.d(2,:,22,:))];
%         [squeeze(samples.d(3,:,14,:)); squeeze(samples.d(3,:,18,:)); squeeze(samples.d(3,:,22,:))]], 1500);

mch      = [1*ones(nt,1), squeeze([samples.d(1,:,14,:), samples.d(2,:,14,:), samples.d(3,:,14,:)]); % P1
            2*ones(nt,1), squeeze([samples.d(1,:,18,:), samples.d(2,:,18,:), samples.d(3,:,18,:)]); % P2
            3*ones(nt,1), squeeze([samples.d(1,:,22,:), samples.d(2,:,22,:), samples.d(3,:,22,:)])]; % P3
post_dmch = mch(:,2:end); %datasample(mch(:,2:end), 1500);


    % Misaligned_incongruent_low        
%     post_d23 = datasample([[squeeze(samples.d(1,:,15,:)); squeeze(samples.d(1,:,19,:)); squeeze(samples.d(1,:,23,:))];
%         [squeeze(samples.d(2,:,15,:)); squeeze(samples.d(2,:,19,:)); squeeze(samples.d(2,:,23,:))];
%         [squeeze(samples.d(3,:,15,:)); squeeze(samples.d(3,:,19,:)); squeeze(samples.d(3,:,23,:))]], 1500);
    
mil      = [1*ones(nt,1), squeeze([samples.d(1,:,15,:), samples.d(2,:,15,:), samples.d(3,:,15,:)]); % P1
            2*ones(nt,1), squeeze([samples.d(1,:,19,:), samples.d(2,:,19,:), samples.d(3,:,19,:)]); % P2
            3*ones(nt,1), squeeze([samples.d(1,:,23,:), samples.d(2,:,23,:), samples.d(3,:,23,:)])]; % P3
post_dmil = mil(:,2:end); %datasample(mil(:,2:end), 1500);

    % Misaligned_incongruent_high
%     post_d24 = datasample([[squeeze(samples.d(1,:,16,:)); squeeze(samples.d(1,:,20,:)); squeeze(samples.d(1,:,24,:))];
%         [squeeze(samples.d(2,:,16,:)); squeeze(samples.d(2,:,20,:)); squeeze(samples.d(2,:,24,:))];
%         [squeeze(samples.d(3,:,16,:)); squeeze(samples.d(3,:,20,:)); squeeze(samples.d(3,:,24,:))]], 1500);

mih      = [1*ones(nt,1), squeeze([samples.d(1,:,16,:), samples.d(2,:,16,:), samples.d(3,:,16,:)]); % P1
            2*ones(nt,1), squeeze([samples.d(1,:,20,:), samples.d(2,:,20,:), samples.d(3,:,20,:)]); % P2
            3*ones(nt,1), squeeze([samples.d(1,:,24,:), samples.d(2,:,24,:), samples.d(3,:,24,:)])]; % P3
post_dmih = mih(:,2:end); %datasample(mih(:,2:end), 1500);

%%
                               
subplot(2,2,1)
fh11 = fill([xt, fliplr(xt)], [low11, fliplr(hi11)], [.27 .47 .77], 'FaceAlpha', .4, 'EdgeAlpha', 0); hold on
fh12 = fill([xt, fliplr(xt)], [low12, fliplr(hi12)], [.56 .03 .52], 'FaceAlpha', .4, 'EdgeAlpha', 0);
fh13 = fill([xt, fliplr(xt)], [low13, fliplr(hi13)], [.10 .72 .88], 'FaceAlpha', .4, 'EdgeAlpha', 0);
fh14 = fill([xt, fliplr(xt)], [low14, fliplr(hi14)], [.43 .27 .66], 'FaceAlpha', .4, 'EdgeAlpha', 0);
vh11 = violin(post_dacl, 'x', xloc, ...
    'facecolor', [.27 .47 .77], 'edgecolor', 'k', 'facealpha', .7); hold on %C-L
vh12 = violin(post_dach, 'x', xloc, ...
    'facecolor', [.56 .03 .52], 'edgecolor', 'k', 'facealpha', .7); %C-H
vh13 = violin(post_dail, 'x', xloc, ...
    'facecolor', [.10 .72 .88], 'edgecolor', 'k', 'facealpha', .7); %I-L
vh14 = violin(post_daih, 'x', xloc, ...
    'facecolor', [.43 .27 .66], 'edgecolor', 'k', 'facealpha', .7); %I-H

set(gca, 'XLim', [0 2.5], 'YLim', [0 5], 'XTick', 0:.5:2.5, 'XTickLabel',  {'.00', '.50', '1.0', '1.5', '2.0', '2.5'})
% set(gca, 'XLim', [0 xmax], 'YLim', [0 5], 'TickLabelInterpreter', 'tex', 'XTick', xloc, 'XTickLabel',  xticklabs) % num2str(round(tpt,2)')
% set(gca, 'Children', [vh11, vh12, fh11, fh12])
xlabel('Total Processing Time (sec)')
ylabel('d-prime')
legend([fh12, fh14, fh11, fh13],...
    'Aligned: Congruent High','Aligned: Incongruent High', 'Aligned: Congruent Low','Aligned: Incongruent Low',...
    'Location', 'Best')
% .44 .99 .69; .54 .76 .60; .67 .53 .54; .87 .25 .45
subplot(2,2,2)
fh21 = fill([xt, fliplr(xt)], [low21, fliplr(hi21)], [.54 .76 .60], 'FaceAlpha', .5, 'EdgeAlpha', 0); hold on
fh22 = fill([xt, fliplr(xt)], [low22, fliplr(hi22)], [.87 .25 .45], 'FaceAlpha', .5, 'EdgeAlpha', 0);
fh23 = fill([xt, fliplr(xt)], [low23, fliplr(hi23)], [.44 .99 .69], 'FaceAlpha', .5, 'EdgeAlpha', 0);
fh24 = fill([xt, fliplr(xt)], [low24, fliplr(hi24)], [.67 .53 .54], 'FaceAlpha', .5, 'EdgeAlpha', 0);

vh21 = violin(post_dmcl, 'x', xloc, ...
    'facecolor', [.54 .76 .60], 'edgecolor', 'k', 'facealpha', .5); hold on
vh22 = violin(post_dmch, 'x', xloc, ...
    'facecolor', [.87 .25 .45], 'edgecolor', 'k', 'facealpha', .5);
vh23 = violin(post_dmil, 'x', xloc, ...
    'facecolor', [.44 .99 .69], 'edgecolor', 'k', 'facealpha', .5);
vh24 = violin(post_dmih, 'x', xloc, ...
    'facecolor', [.67 .53 .54], 'edgecolor', 'k', 'facealpha', .5);

set(gca, 'XLim', [0 2.5], 'YLim', [0 5], 'XTick', 0:.5:2.5, 'XTickLabel',  {'.00', '.50', '1.0', '1.5', '2.0', '2.5'})
% set(gca, 'XLim', [0 xmax], 'YLim', [0 5],'TickLabelInterpreter', 'tex', 'XTick', xloc, 'XTickLabel',  xticklabs)
% set(gca, 'Children', [vh21, vh22, fh21, fh22])
xlabel('Total Processing Time (sec)')
ylabel('d-prime')
legend([fh22, fh24, fh21, fh23],...
    'Misaligned: Congruent High', 'Misaligned: Incongruent High', 'Misaligned: Congruent Low', 'Misaligned: Incongruent Low',...
    'Location', 'Best')

%% Compute overlap between m distributions
sets = allcomb(1:size(mf,1), 1:size(mf,1));
for i = 1:size(sets, 1) 
    olap(sets(i,1), sets(i,2)) = overlap(mf(sets(i,1), :), mf(sets(i,2), :), xi); 
end

%% Compute interaction posterior
% Low salience
figure('WindowStyle', 'docked')
subplot(2,1,1)
lowI = (post_dacl - post_dail) - (post_dmcl - post_dmil);
violin(lowI, 'x', xloc, ...
    'facecolor', [.97 .93 .22], 'edgecolor', 'k', 'facealpha', .7); hold on
plot(xt, zeros(numel(xt),1), '--k')
set(gca, 'XLim', [0 2.5], 'YLim', [-4 4], 'XTick', 0:.5:2.5, 'XTickLabel',  {'.00', '.50', '1.0', '1.5', '2.0', '2.5'})
%set(gca, 'XLim', [0 xmax], 'YLim', [-5 5],'TickLabelInterpreter', 'tex', 'XTick', xloc, 'XTickLabel',  xticklabs)
xlabel('Total Processing Time (sec)')
ylabel('Congruency x Alignment')
%ylabel('Congruency x Alignment: (AC - AI) - (MC - MI)')
title('Low Discriminability')

fprintf('95%s HDI''s for Low Salience Interaction\n', '%')
[prctile(lowI, 2.5); prctile(lowI, 97.5)]

subplot(2,1,2)
hiI = (post_dach - post_daih) - (post_dmch - post_dmih);
violin(hiI, 'x', xloc, ...
    'facecolor', [.96 .73 .27], 'edgecolor', 'k', 'facealpha', .7); hold on
plot(xt, zeros(numel(xt),1), '--k')
set(gca, 'XLim', [0 2.5], 'YLim', [-4 4], 'XTick', 0:.5:2.5, 'XTickLabel',  {'.00', '.50', '1.0', '1.5', '2.0', '2.5'})
%set(gca, 'XLim', [0 xmax], 'YLim', [-5 5],'TickLabelInterpreter', 'tex', 'XTick', xloc, 'XTickLabel',  xticklabs)
xlabel('Total Processing Time (sec)')
ylabel('Congruency x Alignment')
%ylabel('Congruency x Alignment: (AC - AI) - (MC - MI)')
title('High Discriminability')

fprintf('95%s HDI''s for High Salience Interaction\n', '%')
[prctile(hiI, 2.5); prctile(hiI, 97.5)]

%% Parameter interactions
mlow = (mset(:,1) - mset(:,3)) - (mset(:,5) - mset(:,7));
mhi = (mset(:,2) - mset(:,4)) - (mset(:,6) - mset(:,8));
T0low = (T0set(:,1) - T0set(:,3)) - (T0set(:,5) - T0set(:,7));
T0hi = (T0set(:,2) - T0set(:,4)) - (T0set(:,6) - T0set(:,8));
taulow = (tauset(:,1) - tauset(:,3)) - (tauset(:,5) - tauset(:,7));
tauhi = (tauset(:,2) - tauset(:,4)) - (tauset(:,6) - tauset(:,8));

figure('WindowStyle', 'docked')

% mparm
subplot(1,3,1)
xi = linspace(-3, 3, 200);
mLowDens = ksdensity(mlow, xi); 
mHiDens = ksdensity(mhi, xi);
fill([xi, fliplr(xi)], [mLowDens, zeros(1,size(mLowDens,2))], colours(1,:), 'FaceAlpha', .7); hold on
fill([xi, fliplr(xi)], [mHiDens, zeros(1,size(mHiDens,2))], colours(2,:), 'FaceAlpha', .7); hold on
xlabel('m Interaction')
ylabel('Posterior Density')
legend('m Low', 'm High', 'Location', 'Best')

mLow95 = [prctile(mlow, 2.5), prctile(mlow, 97.5)];
fprintf('95%s HDI, m parm, low salience = (%3.2f, %3.2f)\n', '%', mLow95(1), mLow95(2))
mHi95 = [prctile(mhi, 2.5), prctile(mhi, 97.5)];
fprintf('95%s HDI, m parm, high salience = (%3.2f, %3.2f)\n\n', '%', mHi95(1), mHi95(2))

% T0 parm
subplot(1,3,2)
xi = linspace(-3, 3, 200);
T0LowDens = ksdensity(T0low, xi); 
T0HiDens = ksdensity(T0hi, xi);
fill([xi, fliplr(xi)], [T0LowDens, zeros(1,size(T0LowDens,2))], colours(1,:), 'FaceAlpha', .7); hold on
fill([xi, fliplr(xi)], [T0HiDens, zeros(1,size(T0HiDens,2))], colours(2,:), 'FaceAlpha', .7); hold on
xlabel('T0 Interaction')
ylabel('Posterior Density')
legend('T0 Low', 'T0 High', 'Location', 'Best')

T0Low95 = [prctile(T0low, 2.5), prctile(T0low, 97.5)];
fprintf('95%s HDI, T0 parm, low salience = (%3.2f, %3.2f)\n', '%', T0Low95(1), T0Low95(2))
T0Hi95 = [prctile(T0hi, 2.5), prctile(T0hi, 97.5)];
fprintf('95%s HDI, T0 parm, high salience = (%3.2f, %3.2f)\n\n', '%', T0Hi95(1), T0Hi95(2))


% tau parm
subplot(1,3,3)
xi = linspace(-3, 3, 200);
tauLowDens = ksdensity(taulow, xi); 
tauHiDens = ksdensity(tauhi, xi);
fill([xi, fliplr(xi)], [tauLowDens, zeros(1,size(tauLowDens,2))], colours(1,:), 'FaceAlpha', .7); hold on
fill([xi, fliplr(xi)], [tauHiDens, zeros(1,size(tauHiDens,2))], colours(2,:), 'FaceAlpha', .7); hold on
xlabel('tau Interaction')
ylabel('Posterior Density')
legend('tau Low', 'tau High', 'Location', 'Best')

tauLow95 = [prctile(taulow, 2.5), prctile(taulow, 97.5)];
fprintf('95%s HDI, tau parm, low salience = (%3.2f, %3.2f)\n', '%', tauLow95(1), tauLow95(2))
tauHi95 = [prctile(tauhi, 2.5), prctile(tauhi, 97.5)];
fprintf('95%s HDI, tau parm, high salience = (%3.2f, %3.2f)\n\n', '%', tauHi95(1), tauHi95(2))


%% Individual subject interaction posteriors
sub_sets = {[1], [2], [3], [1:3]};
subjectIDs = {'    1', '    2', '    3', 'Group'};
cnt = 1; 
for i = 1:numel(sub_sets)
   lowI = (acl(ismember(acl(:,1), sub_sets{i}), 2:end) -...
              ail(ismember(ail(:,1), sub_sets{i}), 2:end)) -...
          (mcl(ismember(mcl(:,1), sub_sets{i}),2:end) -...
              mil(ismember(mil(:,1), sub_sets{i}), 2:end));

   highI = (ach(ismember(ach(:,1), sub_sets{i}), 2:end) -...
              aih(ismember(aih(:,1), sub_sets{i}), 2:end)) -...
           (mch(ismember(mch(:,1), sub_sets{i}),2:end) -...
              mih(ismember(mih(:,1), sub_sets{i}), 2:end));          

   stats.mean.lowI(i,:) = mean(lowI);       
   stats.mean.highI(i,:) = mean(highI);

   stats.ci_low.lowI(i,:) = prctile(lowI, 2.5);
   stats.ci_low.highI(i,:) = prctile(highI, 2.5);

   stats.ci_high.lowI(i,:) = prctile(lowI, 97.5);
   stats.ci_high.highI(i,:) = prctile(highI, 97.5);

   % Make table
    itemConditionIds = {' Low', 'High'};
    totalProcessingTimeStr = cellstr(num2str(round(tpt * 1000, 0)'))';
    cons = [1 1; 1 2];
    means = [stats.mean.lowI(i,:); stats.mean.highI(i,:)];
    ci_lo = [stats.ci_low.lowI(i,:); stats.ci_low.highI(i,:)];
    ci_hi = [stats.ci_high.lowI(i,:); stats.ci_high.highI(i,:)];

    for j = 1:2
        itable{cnt,1} = sprintf('%s\t %s\t', subjectIDs{i}, itemConditionIds{j});
        for k = 1:D
            itable{cnt,1} = [itable{cnt,1}, sprintf('%1.2f [%1.2f-%1.2f]\t', means(j,k), ci_lo(j,k), ci_hi(j,k))];
        end
        cnt = cnt + 1; 
    end
end





