%% analyse weighted model outputs
% written by danyal akarca, university of cambridge, 2022

% clear
clear; clc;

% add paths
addpath('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/models/matlab');
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer');
addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps');

% load empirical data
load '/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/prepare/data/consensus.mat' consensus;

% change directory
directory = '/imaging/astle/users/da04/Postdoc/weighted_gm/model_outputs_141222';
cd(directory);

%% load all networks

% get list
list = dir(fullfile(directory,'*.mat'));
list = string({list.name});
nnet = length(list);
% initialise key outputs
alpha = zeros(nnet,1);
run = zeros(nnet,1);
energy = zeros(nnet,2);
nks = 4;
ks = zeros(nnet,nks,2);
ks_labels = string({'Strength','Clustering','Betweenness','Edge length'});
% load all networks
for net = 1:nnet;
    % load the network
    load(list(net));
    % extract model related variables
    eta(net) = output.model.settings.eta;
    gamma(net) = output.model.settings.gamma;
	omega(net) = output.model.optimisation.omega;
    alpha(net) = output.model.optimisation.alpha;
    run(net) = output.run;
    % extract evaluations
    energy(net,:) = [output.evaluation.binary.energy, output.evaluation.weighted.energy];
    ks(net,:,1) = [output.evaluation.binary.ks]; ks(net,:,2) = [output.evaluation.weighted.ks];
    % clear output
    clear output
    % display
    disp(sprintf('Network %g of %g loaded',net,nnet));
end
% place all parameters together so that they can be indexed
parameter_settings = [eta' gamma' omega' alpha];
% get indices for networks over the same omega, alpha and run
[~,vo,io] = unique(omega);
[~,va,ia] = unique(alpha);
[~,vr,ir] = unique(run+1);
% number of unique alpha and runs
nomega = length(vo);
nalpha = length(va);
nrun = length(vr);

%% plot energy over the whole sample
% visualise
h = figure; h.Position = [100 100 1600 180];
subplot(1,7,1);
histogram(energy(:,1),'edgecolor','w'); xlim([0 1]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 12; box off;
xlabel('Energy'); ylabel('Frequency');
subplot(1,7,2);
histogram(energy(:,2),'edgecolor','w'); xlim([0 1]);
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 12; box off;
xlabel('Energy'); ylabel('Frequency');
subplot(1,7,3);
scatter(energy(:,1),energy(:,2),50,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 12; box off;
xlabel('Binary energy'); ylabel('Weighted energy'); axis equal;
subplot(1,7,4);
scatter(ks(:,1,1),ks(:,1,2),50,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 12; box off;
xlabel('Binary KSk'); ylabel('Weighted KSk'); axis equal;
subplot(1,7,5);
scatter(ks(:,2,1),ks(:,2,2),50,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 12; box off;
xlabel('Binary KSc'); ylabel('Weighted KSc'); axis equal;
subplot(1,7,6);
scatter(ks(:,3,1),ks(:,3,2),50,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 12; box off;
xlabel('Binary KSb'); ylabel('Weighted KSb'); axis equal;
subplot(1,7,7);
scatter(ks(:,4,1),ks(:,4,2),50,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 12; box off;
xlabel('Binary KSd'); ylabel('Weighted KSd'); axis equal;

%% plot the energy and ks landscapes

% average the outcome at all unique combinations of eta, gamma and omega, alpha
ue = unique(eta); ug = unique(gamma); uo = unique(omega); ua = unique(alpha');
eta_gam = combvec(ue,ug);
ome_alp = combvec(uo,ua);

% find where these are, and average within indices
energy_av_eta_gam = [];
energy_av_ome_alp = [];
for i = 1:length(eta_gam);
    ind = find(sum([eta' gamma']==eta_gam(:,i)',2)==2);
    for bw = 1:2;
        energy_av_eta_gam(i,bw) = mean(energy(ind,bw));
    end
end
for i = 1:length(ome_alp);
    ind = find(sum([omega' alpha]==ome_alp(:,i)',2)==2);
    for bw = 1:2; % loop over binary and weighted energy
        energy_av_ome_alp(i,bw) = mean(energy(ind,bw));
    end
end

% currently it plots a single instance
% best would do it by reshaping the energy matrix!

% plot the weighted landsscape
h = figure; h.Position = [100 100 800 300];
subplot(1,2,1);
scatter(eta_gam(1,:)',eta_gam(2,:)',5000,energy_av_eta_gam(:,1),'.'); c = colorbar; c.Label.String = 'Binary energy';;
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('\eta'); ylabel('\gamma'); caxis([0 1]); colormap(parula); xlim([-3 -2.2]); ylim([0.35 0.4]); % to update
subplot(1,2,2);
scatter(ome_alp(1,:)',ome_alp(2,:)',5000,energy_av_ome_alp(:,1),'.'); c = colorbar; c.Label.String = 'Binary energy';
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('\omega'); ylabel('\alpha'); caxis([0 1]); colormap(parula); xlim([0.8 1.01]); ylim([0.02 0.1]);

% plot the weighted landscape
h = figure; h.Position = [100 100 800 300];
subplot(1,2,1);
scatter(eta_gam(1,:)',eta_gam(2,:)',5000,energy_av_eta_gam(:,2),'.'); c = colorbar; c.Label.String = 'Weighted energy';;
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('\eta'); ylabel('\gamma'); caxis([0 1]); colormap(parula); xlim([-3 -2.2]); ylim([0.35 0.4]); % to update
subplot(1,2,2);
scatter(ome_alp(1,:)',ome_alp(2,:)',5000,energy_av_ome_alp(:,2),'.'); c = colorbar; c.Label.String = 'Weighted energy';
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('\omega'); ylabel('\alpha'); caxis([0 1]); colormap(parula); xlim([0.8 1.01]); ylim([0.02 0.1]);

%% compare distributions to the empirical

% plot an example empirical comparison
% find the network with the lowest energy
[i v] = min(energy);
% select a network
networkv1 = v(2);
% note, or pick another network with different parameters
networkv2 = find(parameter_settings(:,1)==ue(1)&...
    parameter_settings(:,2)==ug(1)&...
    parameter_settings(:,3)==uo(1)&...
    parameter_settings(:,4)==ua(1));

% set 
network = networkv1;

% load final weighted simulation
load(list(network));
Wsynth = squeeze(output.network.weighted(:,:,end));
% get its settings

% load the empirical
Wtgt = squeeze(consensus.dk.connectivity(1,6,:,:));

% plot the empirical statistics
x = cell(4,1);
x{1} = sum(weight_conversion(Wtgt,'normalize'),2);
x{2} = real(clustering_coef_wu(weight_conversion(Wtgt,'normalize')));
x{3} = betweenness_wei(weight_conversion(Wtgt,'normalize'));
D = consensus.dk.euclidean;
x{4} = D(triu(Wtgt,1) > 0);

y = cell(4,1);
y{1} = sum(weight_conversion(Wsynth,'normalize'),2);
y{2} = real(clustering_coef_wu(weight_conversion(Wsynth,'normalize')));
y{3} = betweenness_wei(weight_conversion(Wsynth,'normalize'));
D = consensus.dk.euclidean;
y{4} = D(triu(Wsynth,1) > 0);

% do cdf plots
h = figure; h.Position = [100 100 1200 250];
for stat = 1:4;
    % take the statistics
    x1 = x{stat};
    x2 = y{stat};
    % plot a cfds
    % take common bins
    binEdges =  [-inf;sort([x1;x2]);inf];
    % compute the bin counts
    binCounts1 = histc (x1,binEdges,1);
    binCounts2 = histc (x2,binEdges,1);
    % compute a cumulative probability
    sumCounts1 = cumsum(binCounts1)./sum(binCounts1);
    sumCounts2 = cumsum(binCounts2)./sum(binCounts2);
    % compute the cumulative density
    sampleCDF1 = sumCounts1(1:end-1);
    sampleCDF2 = sumCounts2(1:end-1);
    % compute the diffference 
    deltaCDF = abs(sampleCDF1-sampleCDF2);
    % compute the largest difference (KS)
    kstat = max(deltaCDF);
    % keep
    kstat_instance(stat) = kstat;
    % visualise the CDFs
    subplot(1,4,stat);
    plot(sampleCDF1,'linewidth',2);
    hold on;
    plot(sampleCDF2,'linewidth',2);
    box off;
    xlabel(ks_labels{stat}); ylabel(sprintf('F(%s)',ks_labels(stat)));
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
    title(sprintf('KS=%.3g',kstat_instance(stat)),'FontName','Arial','FontWeight','normal','FontSize',14);
end
sgtitle(sprintf('Distributions. Simulation eta=%.3g, gamma=%.3g, omega=%.3g alpha=%.3g',...
    output.model.settings.eta,output.model.settings.gamma,output.model.optimisation.omega, output.model.optimisation.alpha))

% labels
ks_labels = string({'s','wc','wb','d'});
% visualise
h = figure; h.Position = [10 10 1200 250];
for i = 1:4;
    subplot(1,4,i); 
    histogram(x{i},'edgecolor','w');
    hold on;
    histogram(y{i},'edgecolor','w');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
    xlabel(ks_labels(i)); ylabel('Frequency');
    title(sprintf('KS=%.3g',kstat_instance(i)),'FontName','Arial','FontWeight','normal','FontSize',14);
end
sgtitle(sprintf('Distributions. Simulation eta=%.3g, gamma=%.3g, omega=%.3g alpha=%.3g',...
    output.model.settings.eta,output.model.settings.gamma,output.model.optimisation.omega, output.model.optimisation.alpha))

% visualise the developmental trajectory of this network (e.g., PCA)
trajectory = output.network.weighted;
nsteps = size(trajectory,3);
nstat = 4;
step_stats = zeros(nsteps,nstat);
stat_labels = {'Strength','Efficiency','Modularity','\sigma'};
nperm = 100;
for step = 1:size(trajectory,3);
    net = squeeze(trajectory(:,:,step));
    % compute the weight
    step_stats(step,1) = sum(net,'all');
        % compute the efficiency
    step_stats(step,2) = efficiency_wei(net);
    % compute the modularity
    [~,step_stats(step,3)] = modularity_und(net);
    % compute the small worldness
        % binary version
        A = double(net>0);
        % compute nnode
        nnode = size(A,1);
        % compute number of edges
        m = nnz(A)/2;
        % compute observed
        clu = mean(clustering_coef_bu(A));
        cpl = charpath(A);
        % initialise
        clu_perm = [];
        cpl_perm = [];
        % compute nulls
        for j = 1:nperm;
            % form a random
            Arand = makerandCIJ_und(nnode,m);
            % form a lattic
            clu_perm(j) = mean(clustering_coef_bu(Arand));
            cpl_perm(j) = charpath(Arand);
        end
        % calclate means
        mean_clu_perm = mean(clu_perm);
        mean_cpl_perm = mean(cpl_perm);
        % calculate smw
        step_stats(step,4) = [clu/mean_clu_perm]/[cpl/mean_cpl_perm];
        % display
        disp(sprintf('Step %g of %g computed',step,nsteps));
end
% compute latent space if we want to look at it
[coeff score latent tsq exp] = pca(step_stats,'numcomponents',nstat);
% plot the trajectory across the statisics
h = figure; h.Position = [100 100 800 250];
for stat = 1:nstat;
    subplot(1,nstat,stat);
    scatter(1:nsteps,step_stats(:,stat),500,1:nsteps,'.'); colormap(magma);
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
    xlabel('Simulated time'); ylabel(stat_labels{stat});
end

% evaluate the correlation between distance and edge strength
h = figure; h.Position = [100 100 500 250];
[r1 p1] = corr(Wtgt(:),D(:));
[r2 p2] = corr(Wsynth(:),D(:));
scatter(Wtgt(:),D(:),200,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
hold on;
scatter(Wsynth(:),D(:),200,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
legend({sprintf('Observation, r=%.3g, p=%.3g',r1,p1),sprintf('Simulation, r=%.3g, p=%.3g',r2,p2)},'box','off');
xlabel('Euclidean distance'); ylabel('Edge weight');

% evaluate the correlation between degree and strength
h = figure; h.Position = [100 100 500 250];
[r1 p1] = corr(degrees_und(Wtgt)',strengths_und(Wtgt)');
[r2 p2] = corr(degrees_und(Wsynth)',strengths_und(Wsynth)');
scatter(degrees_und(Wtgt),strengths_und(Wtgt),200,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
hold on;
scatter(degrees_und(Wsynth),strengths_und(Wsynth),200,'.');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
legend({sprintf('Observation, r=%.3g, p=%.3g',r1,p1),sprintf('Simulation, r=%.3g, p=%.3g',r2,p2)},'box','off');
xlabel('Number of connections'); ylabel('Node strength');

%{
%% evaluate hub nodes and their feeders

% set network to plot
Aset = Wsynth;

% maximum level
kmax = max(degrees_und(Aset));
% empirical
[r nk ek] = rich_club_bu(Aset,kmax);
% random permutation after 50 rewires
nperm = 1000;
% initalise
r_perm = zeros(nperm,kmax);
nk_perm = zeros(nperm,kmax);
ek_perm = zeros(nperm,kmax);
for perm = 1:nperm;
    % permute
    Aconses_bin_rand = randmio_und(Aset,50);
    % calculate the permutated rich club
    [r_perm(perm,:) nk_perm(perm,:) ek_perm(perm,:)] = ...
        rich_club_bu(Aconses_bin_rand,kmax);
    % display
    disp(sprintf('Permutation %g of %g completed',perm,nperm));
end
% take ratio
rratio = r./mean(r_perm,1);
% pvalue
prank = nan(kmax,1);
for lvl = 1:kmax;
    % take observed
    a = r(lvl);
    % take permutation
    b = r_perm(:,lvl);
    % rank
    prank(lvl) = sum(a<b)./nperm;
end

% plot the outcome
h = figure; h.Position = [100 100 1200 300];
colours = magma;
% rich club coefficient
subplot(1,3,1);
plot(r,'linewidth',3,'color',colours(100,:)); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 16; box off;
xlabel('Degree, k'); ylabel('\phi');
xlim([0 length(r)]); xline(22,'linewidth',3,'color',colours(255,:));
subplot(1,3,2);
plot(rratio,'linewidth',3,'color',colours(100,:)); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 16; box off;
xlabel('Degree, k'); ylabel('\phi/<\phi_{rand}>');
xlim([0 length(r)]); xline(22,'linewidth',5,'color',colours(255,:));
yline(1,'linestyle','--','color',[.5 .5 .5],'linewidth',3);
subplot(1,3,3); 
plot(prank,'linewidth',3,'color',colours(100,:)); 
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 16; box off;
xlabel('Degree, k'); ylabel('p_{perm}');
xlim([0 length(r)]); xline(22,'linewidth',5,'color',colours(255,:));
yline(0.05,'linestyle','--','color',[.5 .5 .5],'linewidth',3);
%}

%% compute the the time courses and statistical distributions across weighted parameters

% I've seen different modularity trends depending on the binary parameters,
% so we will want to sweep the space (eta gamma) for interesting trends in alpha, omega effects on trajectories.

% to plot every network is excessive, so sample from the lowest average energy eta and gamma
% set eta and gamma
eta_select = ue(1); gamma_select = ug(1);
% loop over all combinations within this parameter combination, take distributions and time courses
nnets = length(uo)*length(ua);
ntime = 227;
nstat = 2;
stat_labels = {'Strength','Modularity','Efficiency'};
step_stats = zeros(nnets,ntime,nstat);
step = 1;
% note, we can fix some variables if we wish.
for o = 1:length(uo);
    for a = 1:length(ua);
        % get the network
        network = find(parameter_settings(:,1)==eta_select&...
        parameter_settings(:,2)==gamma_select&...
        parameter_settings(:,3)==uo(o)&...
        parameter_settings(:,4)==ua(a));
        % load this network
        load(list(network));
        % take the network
        trajectory = output.network.weighted;
        % compute the time course statistics
        nsteps = size(trajectory,3);
        for t = 1:size(trajectory,3);
            net = squeeze(trajectory(:,:,t));
            step_stats(step,t,1) = sum(net,'all');
            [~,step_stats(step,t,2)] = modularity_und(net);
        end
        % display
        disp(sprintf('Network %g of %g computed',step,nnets));
        % update interation
        step = step + 1;
    end
end

%% visualise the time courses

% plot the time courses
h = figure; h.Position = [100 100 600 250];
for stat = 1:nstat;
    subplot(1,nstat,stat);
    for i = 1:nnets;
        plot(step_stats(i,:,stat));
        xlabel('Simulated time');
        hold on;
        ylabel(stat_labels{stat});
    end
b = gca; b.TickDir = 'out'; b.FontSize = 14; b.FontName = 'Arial';
box off; 
end

% take the early and late scores and see if they correlate: highlights to what extent that statistics are determined early in network development
timeA = 100;
timeB = 200;
time_statistics = squeeze(step_stats(:,[timeA timeB],:));
h = figure; h.Position = [100 100 600 250];
for stat = 1:nstat;
    subplot(1,nstat,stat);
    scatter(time_statistics(:,1,stat),time_statistics(:,2,stat),500,1:size(time_statistics,1),'.'); colormap(magma);
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
    xlabel('Early development'); ylabel('Late development');
    title(stat_labels{stat},'FontWeight','normal');
end
