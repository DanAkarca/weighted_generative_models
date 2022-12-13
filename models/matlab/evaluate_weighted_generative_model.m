function output = evaluate_weighted_generative_model(eta,gamma,omega,alpha,run);
    % addpaths
    cd('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/');
    addpath('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/models/matlab/');
    addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT');
    addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer');
    addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps');
    % load the data
    load '/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/prepare/data/consensus.mat' consensus;
    % set directory
    savedir = '/imaging/astle/users/da04/Postdoc/weighted_gm/model_outputs_131222c';
    % run weighted generative model
    % target
    type = 1;
    threshold = 6;
    Wtgt = squeeze(consensus.dk.connectivity(type,threshold,:,:));
    Wseed = squeeze(consensus.dk.seed(type,threshold,:,:));
    D = consensus.dk.euclidean;
    % compute observed binary statistics
    Atgt = double(Wtgt>0);
    xb = cell(4,1);
    xb{1} = degrees_und(Atgt)';
    xb{2} = clustering_coef_bu(Atgt);
    xb{3} = betweenness_bin(Atgt)';
    xb{4} = D(triu(Atgt,1) > 0);
    % compute observed weighted statistics
    xw = cell(4,1);
    xw{1} = sum(weight_conversion(Wtgt,'normalize'),2);
    xw{2} = clustering_coef_wu(weight_conversion(Wtgt,'normalize'));
    xw{3} = betweenness_wei(weight_conversion(Wtgt,'normalize'));
    xw{4} = D(triu(Wtgt,1) > 0);
    % binary generative model parameters
    nnode = 68;
    % seed = Wseed;
    % seed = double(Wseed>0);
    seed = zeros(nnode);
    m = nnz(Wtgt)/2;
    mseed = nnz(seed)/2;
    modeltype = 'matching';
    modelvar = {'powerlaw','powerlaw'};
    params = [eta gamma];
    epsilon = 1e-6;
    % weighted genreative model parameters
    weighted_model = struct;
    weighted_model.update = 1; % flag for weighted genreative models optimisation
    weighted_model.start = 1; % when to start optimisation (this is on top of the initial seed network connectivity)
    weighted_model.optimisation.function = 2; % communicability, 1. normalised communicability, 2.
    weighted_model.optimisation.resolution = 0.05; % how many weight alternations are sampled before the gradient is taken
    weighted_model.optimisation.samples = 5; % how many samples are taken to infer the gradient
    weighted_model.optimisation.alpha = alpha; % the update coefficient [must be a smaller range]
    weighted_model.optimisation.omega = omega; % the local-global weight update term
	weighted_model.settings.eta = eta;
	weighted_model.settings.gamma = gamma;
    % run the weighted genreative model
    [bb cc dd] = weight_optimised_generative_model(seed,D,m,modeltype,modelvar,params,epsilon,...
        weighted_model);
    % get the final network
    wfinal = dd(:,:,end);
    % do evaluations
    % compute compute binary
    Abin = double(wfinal>0); % this may not be accurate - check that this is correct
    Kb = [];
    yb = cell(4,1);
    yb{1} = degrees_und(Abin)';
    yb{2} = clustering_coef_bu(Abin);
    yb{3} = betweenness_bin(Abin)';
    yb{4} = D(triu(Abin,1) > 0);
    for j = 1:4
        Kb(j) = fcn_ks(xb{j},yb{j});
    end
    Eb = max(Kb,[],2);
    % compute weighted
    Kw = [];
    yw = cell(4,1);
    yw{1} = sum(weight_conversion(wfinal,'normalize'),2);
    yw{2} = clustering_coef_wu(weight_conversion(wfinal,'normalize'));
    yw{3} = betweenness_wei(weight_conversion(wfinal,'normalize'));
    yw{4} = D(triu(wfinal,1) > 0); % can cause an error if there are negative weights!
    for j = 1:4
        Kw(j) = fcn_ks(xw{j},yw{j});
    end
    Ew = max(Kw,[],2);
    % keep the output
    output = struct;
    output.model = weighted_model;
    output.network.binary = cc;
    output.network.weighted = dd;
    output.evaluation.binary.energy = Eb;
    output.evaluation.binary.ks = Kb;
    output.evaluation.weighted.energy = Ew;
    output.evaluation.weighted.ks = Kw;
	ouput.detail = 'Evaluated where both the target and simulated network were normalized simultaneously';
    output.run = run;
    % save the output
    savefile = sprintf('%s/output_%g_%g_%g.mat',savedir,run,alpha,omega);
    save(savefile,'output');
    %% ks function
    function kstat = fcn_ks(x1,x2)
        binEdges    =  [-inf ; sort([x1;x2]) ; inf];
        binCounts1  =  histc (x1 , binEdges, 1);
        binCounts2  =  histc (x2 , binEdges, 1);
        sumCounts1  =  cumsum(binCounts1)./sum(binCounts1);
        sumCounts2  =  cumsum(binCounts2)./sum(binCounts2);
        sampleCDF1  =  sumCounts1(1:end-1);
        sampleCDF2  =  sumCounts2(1:end-1);
        deltaCDF  =  abs(sampleCDF1 - sampleCDF2);
        kstat = max(deltaCDF);
    end
end
