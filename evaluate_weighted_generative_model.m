function output = evaluate_weighted_generative_model(Wtgt,seed,D,modeltype,modelvar,params,epsilon,weighted_model);
%   EVALUATE_WEIGHTED_GENERATIVE_MODEL Run weighted generative network model code
%
%   Written by Danyal Akarca, University of Cambridge, 2023
%   For: A weighted generative model of the human connectome
%
%   [output] = evaluate_weighted_generative_model(A,D,m,modeltype,modelvar,params,weightupdate)
%
%   Generates weighted simulated networks using the model described in the study by Akarca et al (2023) in bioRxiv
%   This code was extended from the previous generative_model() function outlined in Betzel et al (2016) in Neuroimage
%
%   Inputs:
%
%           Wtgt,               weighted target network | (nnode x nnode matrix)
%
%           seed,               binary/weighted network of seed connections | (nnode x nnode matrix)
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
%           output              struct of evaluation outputs
%                .model                                 returns back the weighted_model struct
%                .network.binary                        nnode x nnode x m-mseed x nparam matrix of binary networks
%                .network.weighted                      nnode x nnode x m-mseed x nparam matrix of weighted networks
%                .evaluation.binary.energy              returns the binary energy (nparam x 1) over all tested parameters
%                .output.evaluation.binary.ks           returns the binary ks statistics (nparam x 4) over all tested parameters
%                .output.evaluation.weighted.energy     returns the weighted energy (nparam x 1) over all tested parameters
%                .output.evaluation.weighted.ks         returns the binary energy (nparam x 4) over all tested parameters
%
%   See demo.m code for example usage
%
%   Reference:    
%   Akarca et al (2023) bioRxiv
%   Betzel et al (2016) Neuroimage 124:1054-64.

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
nnode = size(Wtgt,1);
m = nnz(Wtgt)/2; 

% run the weighted genreative model
[bb cc dd] = weighted_generative_model(seed,D,m,modeltype,modelvar,params,epsilon,weighted_model);

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