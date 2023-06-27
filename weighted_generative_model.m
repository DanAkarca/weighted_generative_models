function [bb,cc,dd] = weighted_generative_model(A,D,m,modeltype,modelvar,params,epsilon,weighted_model)

%   WEIGHTED_GENERATIVE_MODEL Run weighted generative network model code
%
%   Written by Danyal Akarca, University of Cambridge, 2023
%   For: A weighted generative model of the human connectome
%
%   [bb,cc,dd] = weighted_generative_model(A,D,m,modeltype,modelvar,params,epsilon,weighted_model)
%
%   Generates weighted simulated networks using the model described in the study by Akarca et al (2023) in bioRxiv
%   This code was extended from the previous generative_model() function outlined in Betzel et al (2016) in Neuroimage
%
%   Inputs:
%
%           A,                  binary/weighted network of seed connections | (nnode x nnode matrix)
%
%           D,                  Euclidean distance/fibre length matrix | (nnode x nnode matrix)
%
%           m,                  number of connections that should be present in final simulated network
%
%           modeltype,          the generative rule 'matching' (new rules to be released) | string
%
%           modelvar,           specifies whether the generative rules are based on power-law or exponential relationship ({'powerlaw'}|{'exponential}) | cell
%
%           params,             matrix of wiring parameters ([eta gamma alpha omega]) at which the model should be evaluated | nparams x 4 matrix
%
%           epsilon,            the baseline probability of forming a particular connection (should be a very small number {default = 1e-5})
%
%           weighted_model      struct of settings (see below) which alter the weight optimisation | struct
%                .update,                        indicator to run the weighted generative model {default = 1}
%                .start,                         how many connections (past the seed) to start weight optimisation {default = 1, at beginning}
%                .optimisation.function,         indicator on what communicability model to run. 1 = communicability. 2 = normalised communicability {default = 2}
%                .optimisation.resolution,       resolution used for computing gradient {default = 0.05}
%                .optimisation.samples,          number of samples taken to infer the gradient (default = 5)
%
%   Outputs:
%
%           bb,                 m x nparam matrix of connections
%
%           cc,                 nnode x nnode x m-mseed x nparam matrix of binary networks
%
%           dd,                 nnode x nnode x m-mseed x nparam matrix of weighted networks
%
%   See demo.m code for example usage
%
%   Reference:    
%   Akarca et al (2023) bioRxiv
%   Betzel et al (2016) Neuroimage 124:1054-64.

% if there is no default epsilon, initialise it
if ~exist('epsilon','var')
    epsilon = 1e-5;
end

% compute the number of nodes
n = length(D);

% compute the number of parameters
nparams = size(params,1);

% compute the number of connections in the seed
mseed = nnz(A)/2;

% initialise the output variables
bb = zeros(m,nparams);
cc = zeros(n,n,m-mseed,nparams);
dd = zeros(n,n,m-mseed,nparams);

% set the wiring model type
switch modeltype
    case 'matching' % currently only matching available
        Kseed = matching_ind(A);
        Kseed = Kseed + Kseed';
        for iparam = 1:nparams
            tic
            % set parameters for this run
            eta = params(iparam,1);
            gam = params(iparam,2);
            alpha = params(iparam,3);
            omega = params(iparam,4);
            % run this generative model
            [bb(:,iparam) cc(:,:,:,iparam) dd(:,:,:,iparam)] = fcn_weighted_matching(A,Kseed,D,m,eta,gam,alpha,omega,modelvar,epsilon,weighted_model);
            time = toc;
        end
end

% compute the weighted matching algorithm
function [bb cc dd] = fcn_weighted_matching(A,K,D,m,eta,gam,alpha,omega,modelvar,epsilon,weighted_model)

% take settings from the struct
weight_update = weighted_model.update; % indicator to run the model
start_thr = weighted_model.start; % when the model starts
func = weighted_model.optimisation.function; % which optimisation function
nu = weighted_model.optimisation.resolution; % gradient parameters
rep = weighted_model.optimisation.samples;

% start the model
K = K + epsilon; % add the epsilon
n = length(D); % take the nnode
mseed = nnz(A)/2; % number of edges in seed
mv1 = modelvar{1}; % take if power law or exponential
mv2 = modelvar{2};
% compute the parameterised costs and values for wiring
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
% compute the initial wiring probability 
Ff = Fd.*Fk.*~A; % for non-extant edges
[u,v] = find(triu(ones(n),1)); % compute indices
indx = (v - 1)*n + u;
P = Ff(indx); % get the probability vector
% initalise network simulation
Akeep = zeros(length(A),length(A),m-mseed);
Wkeep = zeros(length(A),length(A),m-mseed);
% loop over added connections
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
        if ii>=mseed+start_thr; % point at which optimisation starts
            % if first iteration, use the biniary
            if ii==mseed+start_thr;
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
                    % can constain the updates to the giant component
                    if func == 1; % calculate the communicability
                        comm = expm(wsynth);
                    else if func == 2; % calculate the normalised communicability
                        s = sum(wsynth,2); % get the strength
                        s(s==0) = epsilon; % if the diagonals are zero, it will fail so provide an epsilon
                        S = diag(s); % diagonalise
                        adj = (S^-.5)*wsynth*(S^-.5); % compute the exponent
                        comm = expm(adj); % calculate the normalised communicability
                        end
                    end
                    % and the summed communicability
                    sum_comm(edge,ru) = sum(comm,'all');
                end
            end
            % compute the gradient
            curve = [];
            for edge = 1:nedge;
                % get the index of this edge
                ind = [a(edge) b(edge)];
                % compute the gradient
                 x = 1:nreps; % get the incremental changes
                 y = (sum_comm(edge,:) * D(ind(1),ind(2))).^omega;
                 z =  fit(x',y','poly1');
                 curve(edge) = z.p1; % first order
            end
            % update the network partially according to the gradient
            for edge = 1:nedge
                % updated optimisation along the gradient
                uv = W(a(edge),b(edge))-(alpha*curve(edge)); 
                % update network
                W(a(edge),b(edge)) = uv;
                W(b(edge),a(edge)) = uv;
                % prevent negative weights
                W(W<0) = 0;
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
