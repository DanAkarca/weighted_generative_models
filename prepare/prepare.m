%% load, expore and prepare cubric data
% written by danyal akarca, university of cambridge
clear; clc;
%% load required data
load('/imaging/astle/users/da04/Postdoc/weighted_gm/original/derivatives/dk/dk_coordinates');
load('/imaging/astle/users/da04/Postdoc/weighted_gm/original/derivatives/dk/dk_descriptives');
load('/imaging/astle/users/da04/Postdoc/weighted_gm/original/derivatives/de/de_coordinates');
load('/imaging/astle/users/da04/Postdoc/weighted_gm/original/derivatives/de/de_descriptives');
addpath('/imaging/astle/users/da04/PhD/toolboxes/distanceDependent');
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps');
%% set directories
nsub = 88;
demographics = readtable('/imaging/astle/users/da04/Postdoc/weighted_gm/original/original_data/Ckids_age_sex.csv');   
cd('/imaging/astle/users/da04/Postdoc/weighted_gm/original/original_data/for_Danyal/');
subjects = dir; 
subjects = string({subjects.name})';
subjects = subjects(4:end);
%% organise demographic data
sex = string(demographics.sex)=="Male"; sex = sex + 1; % female = 1, male = 2
demographics = [demographics.age sex];
%% organise brain data
directory1 = '/imaging/astle/users/da04/Postdoc/weighted_gm/original/original_data/for_Danyal/';
%{
    This contains:
        - signal fraction, after being submitted to COMMIT and normalised by fibre lengths
        - fibre lengths
%}
directory2 = '/imaging/astle/users/da04/Postdoc/weighted_gm/original/original_data/for_Danyal_2/';
%{
    This contains:
        - streamline matrices, after being submitted to COMMIT
        - signal fraction before being normalized by streamline length
%}
directory3 = '/imaging/astle/users/da04/Postdoc/weighted_gm/original/original_data/for_Danyal_origNOS/';
%{
    This contains:
        - original streamline matrices, before being submitted to COMMIT
%}
% number of data types
ntype = 5;
type_label = string({'fibre length','signal fraction length normalised','streamlines','signal fraction','streamlines unfiltered'});
type_label_short = string({'fl','sf','s','usf','us'});
type_label_original = string({'meanL','sumXbyL_over_meanL','NOS','sumX','original_NOS'});
% nnodes
nnode_dk = 85;
nnode_de = 165;
% initialise
dk_connectivity = zeros(nsub,ntype,nnode_dk,nnode_dk);
de_connectivity = zeros(nsub,ntype,nnode_de,nnode_de);
% loop over subjects and load their data
for sub = 1:nsub;
    cd(directory1);
    dk_connectivity(sub,1,:,:) = table2array(readtable(sprintf('%s%s/DesikanKillany_meanL.csv',directory1,subjects(sub))));
    dk_connectivity(sub,2,:,:) = table2array(readtable(sprintf('%s%s/DesikanKillany_sumXbyL_over_meanL.csv',directory1,subjects(sub))));
    dk_connectivity(sub,3,:,:) = table2array(readtable(sprintf('%s%s/DesikanKillany_NOS.csv',directory2,subjects(sub))));
    dk_connectivity(sub,4,:,:) = table2array(readtable(sprintf('%s%s/DesikanKillany_sumX.csv',directory2,subjects(sub))));
    dk_connectivity(sub,5,:,:) = table2array(readtable(sprintf('%s%s/original_NOS_Desikan_Killany',directory3,subjects(sub))));
    de_connectivity(sub,1,:,:) = table2array(readtable(sprintf('%s%s/Destrieux_meanL.csv',directory1,subjects(sub))));
    de_connectivity(sub,2,:,:) = table2array(readtable(sprintf('%s%s/Destrieux_sumXbyL_over_meanL.csv',directory1,subjects(sub))));
    de_connectivity(sub,3,:,:) = table2array(readtable(sprintf('%s%s/Destrieux_NOS.csv',directory2,subjects(sub))));
    de_connectivity(sub,4,:,:) = table2array(readtable(sprintf('%s%s/Destrieux_sumX.csv',directory2,subjects(sub))));
    de_connectivity(sub,5,:,:) = table2array(readtable(sprintf('%s%s/original_NOS_Destrieux',directory3,subjects(sub))));
    % display
    disp(sprintf('Subject %g of %g data loaded',sub,nsub));
end
%% examine unthresholded data
% loop over subjects and compute the correlation across the data
correlation_subjects_dk = zeros(nsub,ntype,ntype);
correlation_subjects_de = zeros(nsub,ntype,ntype);
for sub = 1:nsub;
    % keep the subjects edge-wise data
    y1 = []; y2 = [];
    for type = 1:ntype;
        x1 = squeeze(dk_connectivity(sub,type,:,:));
        y1(:,type) = x1(:);
        x2 = squeeze(de_connectivity(sub,type,:,:));
        y2(:,type) = x2(:);
    end
    % remove non existant elements
    y1 = y1(any(y1,2),:);
    y2 = y2(any(y2,2),:);
    % correlations
    correlation_subjects_dk(sub,:,:) = corr(y1);
    correlation_subjects_de(sub,:,:) = corr(y2);
end
% visualise the average across data types
h = figure; h.Position = [100 100 1000 350];
subplot(1,2,1);
imagesc(squeeze(mean(correlation_subjects_dk,1)));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
colormap(viridis); c = colorbar; c.Label.String = 'correlation';
title('dk');
xticks(1:ntype); yticks(1:ntype);
xticklabels(type_label_short);  yticklabels(type_label_short);
subplot(1,2,2);
imagesc(squeeze(mean(correlation_subjects_de,1)));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
colormap(viridis); c = colorbar; c.Label.String = 'correlation';
title('de');
xticks(1:ntype); yticks(1:ntype);
xticklabels(type_label_short);  yticklabels(type_label_short);
% plot single subject
subject = 1; % set subject
dataplot = [1 2]; % set x and y axis
h = figure; h.Position = [100 100 1000 350];
subplot(1,2,1);
x = squeeze(dk_connectivity(sub,dataplot(1),:,:));
y = squeeze(dk_connectivity(sub,dataplot(2),:,:));
scatter(x(:),y(:),'.','markeredgecolor','k');
xlabel(type_label(dataplot(1))); ylabel(type_label(dataplot(2)));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
subplot(1,2,2);
x = squeeze(de_connectivity(sub,dataplot(1),:,:));
y = squeeze(de_connectivity(sub,dataplot(2),:,:));
scatter(x(:),y(:),'.','markeredgecolor','k');
xlabel(type_label(dataplot(1))); ylabel(type_label(dataplot(2)));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
%% threshold cortical data
% indices of cortical nodes
dk_cortical = dk_descriptives(:,2)=="Cortical";
de_cortical = de_descriptives(:,2)=="Cortical";
% cortical nnodes
nnode_dk_cortical = 68;
nnode_de_cortical = 148;
% set data to threshold
datathreshold = [2 3 4 5];
ndata = length(datathreshold);
% set enforced average densities 
thresholds = [0.35 0.3 0.25 0.2 0.15 0.1]; % set
nthresholds = length(thresholds);
% set resolution of absolute thresholding (sf and stream thresholded differently)
resolution_sf = linspace(0,0.3,2e2); % signal fraction
resolution_se = linspace(1,900,900); % streamlines
% initialise
dk_conectivity_thresholded = zeros(nsub,ndata,nthresholds,nnode_dk_cortical,nnode_dk_cortical);
dk_threshold_required = zeros(ndata,nthresholds);
dk_density_acheived = zeros(ndata,nthresholds,2);
de_conectivity_thresholded = zeros(nsub,ndata,nthresholds,nnode_de_cortical,nnode_de_cortical);
de_threshold_required = zeros(ndata,nthresholds);
de_density_acheived = zeros(ndata,nthresholds,2);
% dk atlas
% loop over enforced threshold for each subject
for threshold = 1:nthresholds;
    % loop over data type
    for type = 1:ndata;
        % get data
        ti = datathreshold(type);
        % display
        disp(sprintf('Determining threshold to acheive %g density for %s data in the dk atlas...',thresholds(threshold)*100,type_label(ti)));
        % threshold to signal fraction or streamline count depending on data type
        if type == 1 | type == 3
            resolution = resolution_sf;
        else
            resolution = resolution_se;
        end
        % loop over these thresholds
        for res = 1:length(resolution);
            % loop over subjects
            for sub = 1:nsub;
                % get the data
                dk = squeeze(dk_connectivity(sub,ti,:,:));
                % remove self connections
                dk(find(eye(size(dk)))) = 0;
                % keep only cortical
                dk = dk(dk_cortical,dk_cortical);
                % apply threshold 
                d = threshold_absolute(dk,resolution(res));
                density(sub) = density_und(d);
            end
            % compute the mean density
            mean_density = mean(density);
            std_density = std(density);
            % stop thresholding if it hits the threshold
            if mean_density < thresholds(threshold);
                % keep this
                dk_threshold_required(type,threshold) = resolution(res);
                dk_density_acheived(type,threshold,1) = mean_density;
                dk_density_acheived(type,threshold,2) = std_density;
                % then apply the threshold to all subjects
                for sub = 1:nsub;
                    % get the correct data
                    ti = datathreshold(type);
                    % get subject data
                    dk = squeeze(dk_connectivity(sub,ti,:,:));
                    % remove self connections
                    dk(find(eye(size(dk)))) = 0;
                    % only cortical
                    dk = dk(dk_cortical,dk_cortical);
                    % apply
                    dk_connectivity_thresholded(sub,type,threshold,:,:) = threshold_absolute(dk,dk_threshold_required(type,threshold));
                end
            break
            end
        end
    end
end
% de atlas
% loop over enforced threshold for each subject
for threshold = 1:nthresholds;
    % loop over data type
    for type = 1:ndata;
        % get data
        ti = datathreshold(type);
        % display
        disp(sprintf('Determining threshold to acheive %g density for %s data in the de atlas...',thresholds(threshold)*100,type_label(ti)));
        % threshold to signal fraction or streamline count depending on data type
        if type == 1 | type == 3
            resolution = resolution_sf;
        else
            resolution = resolution_se;
        end
        % loop over these thresholds
        for res = 1:length(resolution);
            % loop over subjects
            for sub = 1:nsub;
                % get the data
                de = squeeze(de_connectivity(sub,ti,:,:));
                % remove self connections
                de(find(eye(size(de)))) = 0;
                % keep only cortical
                de = de(de_cortical,de_cortical);
                % apply threshold 
                d = threshold_absolute(de,resolution(res));
                density(sub) = density_und(d);
            end
            % compute the mean density
            mean_density = mean(density);
            std_density = std(density);
            % stop thresholding if it hits the threshold
            if mean_density < thresholds(threshold);
                % keep this
                de_threshold_required(type,threshold) = resolution(res);
                de_density_acheived(type,threshold,1) = mean_density;
                de_density_acheived(type,threshold,2) = std_density;
                % then apply the threshold to all subjects
                for sub = 1:nsub;
                    % get the correct data
                    ti = datathreshold(type);
                    % get subject data
                    de = squeeze(de_connectivity(sub,ti,:,:));
                    % remove self connections
                    de(find(eye(size(de)))) = 0;
                    % only cortical
                    de = de(de_cortical,de_cortical);
                    % apply
                    de_connectivity_thresholded(sub,type,threshold,:,:) = threshold_absolute(de,de_threshold_required(type,threshold));
                end
            break
            end
        end
    end
end
%% explore individual differences in thresholded data
% set data
type = 1;
% set connectivity to compute
connectivity = squeeze(dk_connectivity_thresholded(:,type,:,:,:));
% initialise
nmeasures = 2;
statistics = zeros(nsub,nmeasures);
% compute statistics for each subject
for sub = 1:nsub;
    for threshold = 1:nthresholds;
        % get subject data
        w = squeeze(connectivity(sub,threshold,:,:));
        % compute statistics
        statistics(sub,threshold,1) = sum(w,'all');
        [~,statistics(sub,threshold,2)] = modularity_und(w);
    end
end
% visualise
h = figure; h.Position = [100 100 1200 400];
for i = 1:nmeasures;
    subplot(1,nmeasures,i);
    plot(squeeze(statistics(:,:,i))');
    box off;
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
end
%% generate consensus networks on thresholded data
%{
    generates a group-representative structural connectivity matrix by preserving 
    within-/between-hemisphere connection length distributions
%}
% compute the euclidean
dk_euclidean = squareform(pdist(dk_coordinates));
de_euclidean = squareform(pdist(de_coordinates));
% initialise
dk_connectivity_consensus = zeros(ndata,nthresholds,nnode_dk_cortical,nnode_dk_cortical);
de_connectivity_consensus = zeros(ndata,nthresholds,nnode_de_cortical,nnode_de_cortical);
% dk atlas
% loop over thresholds
for threshold = 1:nthresholds;
    % loop over data type
    for type = 1:ndata;
        % get all subject data
        data = squeeze(dk_connectivity_thresholded(:,type,threshold,:,:));
        % set euclidean distances of node centroids
        % compute nbins as the sqrt of the mean density
        D = [];
        for sub = 1:nsub
            A = squeeze(data(sub,:,:));
            D(sub) = density_und(A);
        end
        nbins = round(sqrt(100.*mean(D)));
        % reorder network for function
        data = permute(data,[2 3 1]);
        % input hemispheric ID
        hemiid = zeros(size(dk_descriptives,1),1);
        hemiid(dk_descriptives(:,3)=="Left" & dk_descriptives(:,2)=="Cortical") = 1;
        hemiid(dk_descriptives(:,3)=="Right" & dk_descriptives(:,2)=="Cortical") = 2;
        hemiid(hemiid==0) = [];
        % compute threhsolding procedure for the binary network
        G = fcn_group_bins(data,dk_euclidean,hemiid,nbins);
        % take the sample meaned cortical aspect of the thresholded weighted network
        Av = squeeze(mean(squeeze(dk_connectivity_thresholded(:,type,threshold,:,:)),1));
        % filter the average by the group consensus
        W = Av .* G;
        % save the data
        dk_connectivity_consensus(type,threshold,:,:) = W;
    end
end
% de atlas
% loop over thresholds
for threshold = 1:nthresholds;
    % loop over data type
    for type = 1:ndata;
        % get all subject data
        data = squeeze(de_connectivity_thresholded(:,type,threshold,:,:));
        % set euclidean distances of node centroids
        % compute nbins as the sqrt of the mean density
        D = [];
        for sub = 1:nsub
            A = squeeze(data(sub,:,:));
            D(sub) = density_und(A);
        end
        nbins = round(sqrt(100.*mean(D)));
        % reorder network for function
        data = permute(data,[2 3 1]);
        % input hemispheric ID
        hemiid = zeros(size(de_descriptives,1),1);
        hemiid(de_descriptives(:,3)=="Left" & de_descriptives(:,2)=="Cortical") = 1;
        hemiid(de_descriptives(:,3)=="Right" & de_descriptives(:,2)=="Cortical") = 2;
        hemiid(hemiid==0) = [];
        % compute threhsolding procedure for the binary network
        G = fcn_group_bins(data,de_euclidean,hemiid,nbins);
        % take the sample meaned cortical aspect of the thresholded weighted network
        Av = squeeze(mean(squeeze(de_connectivity_thresholded(:,type,threshold,:,:)),1));
        % filter the average by the group consensus
        W = Av .* G;
        % save the data
        de_connectivity_consensus(type,threshold,:,:) = W;
    end
end
% visualise the outcomes
% set data
type = 1;
% set dataplot
dataplot = dk_connectivity_consensus;
% initialsie
stat = [];
% visualise
h = figure; h.Position = [100 100 2000 240];
for threshold = 1:nthresholds;
    % get data
    x = squeeze(dataplot(type,threshold,:,:));
    % calculate some global statisics
    stat(1,threshold) = sum(x,'all');
    [~,stat(2,threshold)] = modularity_und(x);
    stat(3,threshold) = sum(expm(x),'all');
    % plot
    subplot(1,nthresholds,threshold);
    imagesc(x);
    xlabel('Node'); ylabel('Node');
    b = gca; b.FontName = 'Arial'; colormap(viridis);
end
%% compute the seed network < to review
% set cortical
dk_connectivity_cortical = dk_connectivity(:,datathreshold,dk_cortical,dk_cortical);
de_connectivity_cortical = de_connectivity(:,datathreshold,de_cortical,de_cortical);
% initialise
dk_connectivity_seed = [];
de_connectivity_seed = [];
for type = 1:ndata;
    % dk atlas
    dk_connectivity_type = squeeze(dk_connectivity_cortical(:,type,:,:));
    % take minimum across subjects
    dk_connectivity_seed(type,:,:) = min(dk_connectivity_type,[],1);
    % de atlas
    de_connectivity_type = squeeze(de_connectivity_cortical(:,type,:,:));
    % take minimum across subjects
    de_connectivity_seed(type,:,:) = min(de_connectivity_type,[],1);
end
%% keep data
consensus = struct;
consensus.de.connectivity = de_connectivity_consensus;
consensus.de.seed = de_connectivity_seed;
consensus.de.coordinates = de_coordinates;
consensus.de.euclidean = de_euclidean;
consensus.de.parcellation = de_descriptives;
consensus.dk.connectivity = dk_connectivity_consensus;
consensus.dk.seed = dk_connectivity_seed;
consensus.dk.coordinates = dk_coordinates;
consensus.dk.euclidean = dk_euclidean;
consensus.dk.parcellation = dk_descriptives;
consensus.info.labels = type_label(datathreshold);
consensus.info.thresholds = thresholds;
consensus.version = string({'written by danyal akarca, 11/22. CUBRIC data n=88'});
save('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/prepare/data/consensus.mat','consensus');

consensus_python = struct;
consensus_python.dk_consensus = squeeze(dk_connectivity_consensus(1,6,:,:));
consensus_python.dk_euclidean = dk_euclidean;
consensus_python.version = string({'written by danyal akarca, 11/22. CUBRIC data n=88'});
save('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/prepare/data/consensus_python.mat','consensus_python');