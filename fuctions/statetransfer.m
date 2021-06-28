function [Gamma, vpath, hmm]=statetransfer(data,T,options,h)

%%

% The script corresponding to the transfer analysis used in manuscript 'Spontaneous transient brain states in EEG source space of disorders of consciousness'
% The underlying details could be found in the supplementary information
% an OSL and HMM toolbox are needed, which can be accessed from https://ohba-analysis.github.io/osl-docs/ and https://github.com/OHBA-analysis/HMM-MAR/wiki
% This script was modified from a pipeline: Vidaurre, D., Hunt, L. T., Quinn, A. J., Hunt, B. a. E., Brookes, M. J., Nobre, A. C. & Woolrich, M. W. 2018. Spontaneous cortical activity transiently organises into frequency specific phase-coupling networks. Nat Commun, 9, 2987.
% The original pipeline on MEG analysis could be found from https://github.com/OHBA-analysis/HMM-MAR/tree/master/examples

% Input: data (ROIs * time points), T (length of each session/epoch of single subject)
% options       structure with the training options - see documentation in 
%                       https://github.com/OHBA-analysis/HMM-MAR/wiki
% OUTPUT:
% hmm           estimated HMMMAR model
% Gamma         Time courses of the states probabilities given data
% vpath         most likely state path of hard assignments


%  edit by Yang Bai 2021-06-22

%%
%%%%%to do the preprocessing includes the PCA of the data
N = length(T);
dat = cell(N,1); TT = cell(N,1);
for i=1:N
  t = 1:T(i);
  dat{i} = data(t,:); TT{i} = T(i);
  try data(t,:) = []; 
  catch, error('The dimension of data does not correspond to T');
  end
end
data = dat; T = TT; clear dat TT
[options,data] = checkoptions(options,data,T,0);
options.As = [];
options.A = highdim_pca(data,T,options.pca,...
options.embeddedlags,options.standardise,...
options.onpower,options.varimax,options.detrend,...
options.filter,options.leakagecorr,options.Fs,options.As);  
options.pca = size(options.A,2);
options.ndim = size(options.A,2);
options.S = ones(options.ndim);
options.Sind = formindexes(options.orders,options.S);
if ~options.zeromean, options.Sind = [true(1,size(options.Sind,2)); options.Sind]; end
options.ndim = size(options.A,2);

[hmm,info] = hmmsinit_S1(data,T,options,h);

[hmm,fehist,feterms,rho] = hmmstrain_S(data,T,hmm,info,options);

Gamma = []; Xi = []; vpath = []; residuals = [];

[Gamma,Xi] = hmmdecode(data,T,hmm,0); 

vpath = hmmdecode(data,T,hmm,1); 

    
end

