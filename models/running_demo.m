%% running demo 
% written by danyal akarca
clear; clc;
cd('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/');
addpath('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/');
addpath('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/MIT_code_2011/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/2019_03_03_BCT/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/colorBrewer/');
addpath('/imaging/astle/users/da04/PhD/toolboxes/Colormaps/Colormaps (5)/Colormaps/');
%% load data
% load the data
load('/imaging/astle/users/da04/Postdoc/weighted_gm/weighted_generative_models/prepare/data/consensus.mat');
%% run weighted generative model
% target
type = 1;
threshold = 6;
Wtgt = squeeze(consensus.dk.connectivity(type,threshold,:,:));
Wseed = squeeze(consensus.dk.seed(type,:,:)); Wseed(Wseed<0.5)=0;
D = consensus.dk.euclidean;
% compute observed statistics
x = cell(4,1);
x{1} = sum(weight_conversion(Wtgt,'normalize'),2);
x{2} = clustering_coef_wu(weight_conversion(Wtgt,'normalize'));
x{3} = betweenness_wei(weight_conversion(Wtgt,'normalize'));
x{4} = D(triu(Wtgt,1) > 0);
% binary generative model parameters
nnode = 68;
A = zeros(nnode,nnode);
m = nnz(Wtgt)/2;
modeltype = 'matching';
modelvar = {'powerlaw','powerlaw'};
params = [-2 0.35];
epsilon = 1e-6;
% weighted genreative model parameters
weighted_model = struct;
weighted_model.update = 1; % flag for weighted genreative models optimisation
weighted_model.omega = 2; % the parameterisation of communicability when optimised
weighted_model.start = 1; % when to start optimisation
weighted_model.optimisation.resolution = 0.05; % how many weight alternations are sampled before the gradient is taken
weighted_model.optimisation.samples = 1; % how many samples are taken to infer the gradient
weighted_model.optimisation.alpha = 1; % the update coefficient
%% run the model open
% initalise
Akeep = zeros(length(A),length(A),m);
Wkeep = zeros(length(A),length(A),m);
% take settings
weight_update = weighted_model.update;
omega = weighted_model.omega;
start_thr = weighted_model.start;
alpha = weighted_model.optimisation.alpha;
nu = weighted_model.optimisation.resolution;
rep = weighted_model.optimisation.samples;
eta = params(1);
gam = params(2);
% initialise
keepcurve = {};
% start simulation
K = matching_ind(A);
K = K + K';
K = K + epsilon;
n = length(D);
mseed = nnz(A)/2;
mv1 = modelvar{1};
mv2 = modelvar{2};
switch mv1
    case 'powerlaw'
        Fd = D.^eta;
    case 'exponential'
        Fd = exp(eta*D);
end
switch mv2
    case 'powerlaw'
        Fk = K.^gam;
    case 'exponential'
        Fk = exp(gam*K);
end
Ff = Fd.*Fk.*~A;
[u,v] = find(triu(ones(n),1));
indx = (v - 1)*n + u;
P = Ff(indx);
c = zeros(n,n,m-(mseed+1));
for ii = (mseed + 1):m
    % display
    disp(sprintf('Connection %g of %g added',ii,m));
    C = [0; cumsum(P)];
    r = sum(rand*C(end) >= C);
    uu = u(r);
    vv = v(r);
    A(uu,vv) = 1;
    A(vv,uu) = 1;
    % keep this update
    Akeep(:,:,ii-mseed) = A;
    % if there is a weight update term
    if weight_update == 1;
        % timing
        t(ii) = tic;
        % start optimisation
        if ii>=start_thr; % point at which optimisation starts
            % if first iteration, use the biniary
            if ii==start_thr;
                W = A;
            else
                W = Wkeep(:,:,ii-mseed-1); % take the previously computed version
                W(uu,vv) = 1; % provide the edge that has been added
                W(vv,uu) = 1; 
            end
            % display
            disp(sprintf('Weights at time step %g of %g tuned',ii,m));
            % calculate the communicability gradient - this will all be parameterised
            repvec = [-rep*nu:nu:0:nu:rep*nu];
            nreps = length(repvec);
            [a b] = find(triu(W));
            nedge = length(a);
            % initialise
            norm_comm = [];
            sum_comm = [];
            % simulate out
            for edge = 1:nedge;
                ind = [a(edge) b(edge)];
                % get its value
                vu = W(ind(1),ind(2));
                % get a vector of delta values
                dv = nu*vu;
                % set the total changes
                reps = [vu-rep*dv:dv:vu:dv:vu+rep*dv];
                % loop over each
                for ru = 1:nreps;
                    % get original weighted network
                    wsynth = W;
                    % update it
                    wsynth(ind(1),ind(2)) = reps(ru);
                    wsynth(ind(2),ind(1)) = reps(ru);
                    % calculate the normalised communicability
                    %{
                    s = sum(wsynth,2); % get the strength
                    s(s==0) = epsilon; % if the diagonals are zero, it will fail so provide an epsilon
                    S = diag(s); % diagonalise
                    adj = (S^-.5)*wsynth*(S^-.5); 
                    comm = expm(adj);
                    %}
                    % calculate the communicability
                    comm = expm(wsynth);
                    % and the summed communicability
                    sum_comm(edge,ru) = sum(comm,'all');
                end
            end
            % compute the gradient
            curve = [];
            for edge = 1:nedge;
                % compute the gradient
                 x = 1:nreps;
                 y = sum_comm(edge,:);
                 z =  fit(x',y','poly1');
                 curve(edge) = z.p1; % first order
            end
            keepcurve{ii} = curve;
            % update the network partially according to the gradient
            for edge = 1:nedge
                % updated optimisation along the gradient
                uv = W(a(edge),b(edge))-alpha*curve(edge); 
                % update network
                W(a(edge),b(edge)) = uv;
                W(b(edge),a(edge)) = uv;
            end
            % keep this network for the next iteration
            Wkeep(:,:,ii-mseed) = W;
        end
    end
    % update the probability distribution for the next iteration
    updateuu = find(A*A(:,uu));
    updateuu(updateuu == uu) = [];
    updateuu(updateuu == vv) = [];
    updatevv = find(A*A(:,vv));
    updatevv(updatevv == uu) = [];
    updatevv(updatevv == vv) = [];
    c1 = [A(:,uu)', A(uu,:)];
    for i = 1:length(updateuu)
        j = updateuu(i);
        c2 = [A(:,j)' A(j,:)];
        use = ~(~c1&~c2);
        use(uu) = 0;  use(uu+n) = 0;
        use(j) = 0;  use(j+n) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            K(uu,j) = epsilon;
            K(j,uu) = epsilon;
        else
            K(uu,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
            K(j,uu) = K(uu,j);
        end
    end
    c1 = [A(:,vv)', A(vv,:)];
    for i = 1:length(updatevv)
        j = updatevv(i);
        c2 = [A(:,j)' A(j,:)];
        use = ~(~c1&~c2);
        use(vv) = 0;  use(vv+n) = 0;
        use(j) = 0;  use(j+n) = 0;
        ncon = sum(c1(use))+sum(c2(use));
        if (ncon==0)
            K(vv,j) = epsilon;
            K(j,vv) = epsilon;
        else
            K(vv,j) = (2*(sum(c1(use)&c2(use))/ncon)) + epsilon;
            K(j,vv) = K(vv,j);
        end
    end
    switch mv2
        case 'powerlaw'
            Fk = K.^gam;
        case 'exponential'
            Fk = exp(gam*K);
    end
    Ff = Fd.*Fk.*~A;
    P = Ff(indx);
end
bb = find(triu(A,1)); % keep the binary matrix
cc = Akeep; % keep the binary matrix over time
dd = Wkeep; % keep the weighted matrix over time
dd(:,:,1:start_thr) = cc(:,:,1:start_thr);
%% look at outputs
% output
wfinal = dd(:,:,end);
figure; imagesc(wfinal);
% come some statistical measures over time
nstat = 2;
statistics = zeros(m,nstat);
for step = 1:m;
    % take network
    w = squeeze(dd(:,:,step));
    % compute statistics
    statistics(step,1) = sum(w,'all'); % total weight
    [~,statistics(step,2)] = modularity_und(w); % modularity
end
% plot the weights over time
h = figure; h.Position = [10 10 800 300];
subplot(1,2,1);
plot(statistics(:,1));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('Simulated time'); ylabel('Network weight'); box off;
subplot(1,2,2);
plot(statistics(:,2));
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
xlabel('Simulated time');  ylabel('Modularity Q'); box off;

% plot the final graph
h = figure; h.Position = [10 10 800 300];
subplot(1,2,1);
imagesc(weight_conversion(Wtgt,'normalize'));
c = colorbar; c.Label.String = 'Weight'; colormap(viridis);
xlabel('Node'); ylabel('Node');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
box off;
subplot(1,2,2);
imagesc(weight_conversion(wfinal,'normalize'));
c = colorbar; c.Label.String = 'Weight'; colormap(viridis);
xlabel('Node'); ylabel('Node');
b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14;
box off;

% plo the empirical statistics
y = cell(4,1);
y{1} = sum(weight_conversion(wfinal,'normalize'),2);
y{2} = real(clustering_coef_wu(weight_conversion(wfinal,'normalize')));
y{3} = betweenness_wei(weight_conversion(wfinal,'normalize'));
y{4} = D(triu(wfinal,1) > 0); % can cause an error if there are negative weights!ot distributions

% labels
ks_labels = string({'s','wc','wb','d'});
% visualise
h = figure; h.Position = [10 10 1500 300];
for i = 1:4;
    subplot(1,4,i); 
    histogram(y{i},10,'edgecolor','w');
    b = gca; b.TickDir = 'out'; b.FontName = 'Arial'; b.FontSize = 14; box off;
    xlabel(ks_labels(i)); ylabel('Frequency');
end

% plot out the changes in curves
figure;
for i = 20:m;
    histogram(keepcurve{i},'facealpha',.1,'edgecolor','w');
    hold on;
end