
function [hmm,info] = hmmsinit_S1(Xin,T,options,h)
%%
% The script was a modification from a HMMMAR toolbox function to initialisation before stochastic HMM variational inference
% Please also refer to the original function: hmminit in the HMMMAR toolbox， provide by Diego Vidaurre, OHBA, University of Oxford (2018)
% INPUTS
% Xin: cell with strings referring to the files containing each subject's data, 
%       or cell with with matrices (time points x channels) with each
%       subject's data
% T: cell of vectors, where each element has the length of each trial per
%       subject. Dimension of T{n} has to be (1 x nTrials)
% options: HMM options for both the subject and the group runs
% h: template hmm model used for state transfer
%
%  edit by Yang Bai 2021-06-22
%%

N = length(T); K = options.K;
X = loadfile(Xin{1},T{1},options); ndim = size(X,2);
subjfe_init = zeros(N,3);
loglik_init = zeros(N,1);
pcaprec = options.pcapred>0;

S = options.S==1; regressed = sum(S,1)>0;
if isfield(options,'B') && ~isempty(options.B)
    npred = length(options.orders)*size(options.B,2) + (~options.zeromean);
elseif pcaprec
    npred = options.pcapred + (~options.zeromean);
else
    npred = length(options.orders)*ndim + (~options.zeromean);
end
info = struct();
if isfield(options,'initial_hmm') 
    initial_hmm = options.initial_hmm; 
    options = rmfield(options,'initial_hmm');
else 
    initial_hmm = [];
end
tp_less = max(options.embeddedlags) + max(-options.embeddedlags);
if options.downsample > 0
    downs_ratio = (options.downsample/options.Fs);
else
    downs_ratio = 1; 
end

% init sufficient statistics
subj_m_init = zeros(npred,ndim,N,K);
subj_gram_init = zeros(npred,npred,N,K);
if strcmp(options.covtype,'diag')
    subj_err_init = zeros(ndim,N,K); 
elseif strcmp(options.covtype,'full')
    subj_err_init = zeros(ndim,ndim,N,K);  
else
    subj_err_init = [];
end
subj_time_init = zeros(N,K);

% init subject parameters
Dir2d_alpha_init = zeros(K,K,N); Dir_alpha_init = zeros(K,N);

best_fe = Inf;
I = randpermNK(N,round(N/options.BIGNinitbatch));
covered = [];
for ii = 1:length(I)
    % read data
    subset = I{ii};
    % subset=1:size(data,1); %%%%%%% subset to define how many trials to be processed together
    if size(subset,1)==1, subset = subset'; end
    covered = [covered; subset];
    [X2,XX,Y,Ti2] = loadfile(Xin(subset),T(subset),options);  %%%%%%%%% in the loadfile
    XX_i = cell(1); XX_i{1} = XX;
    if ii==1
        range_data = range(Y);
    else
        range_data = max(range_data,range(Y));
    end
%%%%%%%%%%%%%----------------------------------------------------- first step, to do a initiation
    [hmm_i,Gamma,Xi] = hmmsinit_S2(X2,Ti2,options,h);
    hmm_i.train.pca = options.pca; 
    hmm_i.train.pca_spatial = options.pca_spatial;
    hmm_i.train.embeddedlags = options.embeddedlags;
    hmm_i.train.filter = options.filter;
    hmm_i.train.detrend = options.detrend;
    hmm_i.train.onpower = options.onpower;
    hmm_i.train.downsample = options.downsample;
    hmm_i.train.useParallel = options.useParallel;
    hmm_i.train.grouping = options.grouping;
%%%%%%%%%%%%%%----------------------------------------------------------------- second step, to train Gamma and state transition,
    if ii==1 % get priors
        Dir2d_alpha_prior = hmm_i.prior.Dir2d_alpha;
        Dir_alpha_prior = hmm_i.prior.Dir_alpha;
        hmm_i.train.orders = formorders(hmm_i.train.order,hmm_i.train.orderoffset,...
            hmm_i.train.timelag,hmm_i.train.exptimelag);
        if options.pcapred>0
            Sind = true(options.pcapred,ndim);
        else
            Sind = formindexes(hmm_i.train.orders,hmm_i.train.S)==1;
        end
        if ~hmm_i.train.zeromean, Sind = [true(1,ndim); Sind]; end
    end
        

    tacc = 0; tacc2 = 0;
    for i=1:length(subset)
        subj = subset(i);
        Tsubj = ceil(downs_ratio*(T{subj}-tp_less)); 
        for trial=1:length(Tsubj)
            t = tacc + 1;
            Dir_alpha_init(:,subj) = Dir_alpha_init(:,subj) + Gamma(t,:)';
            tacc = tacc + Tsubj(trial) - hmm_i.train.maxorder;
            t = tacc2 + (1:Tsubj(trial)-hmm_i.train.maxorder-1);
            Dir2d_alpha_init(:,:,subj) = Dir2d_alpha_init(:,:,subj) + squeeze(sum(Xi(t,:,:),1));
            tacc2 = tacc2 + Tsubj(trial) - hmm_i.train.maxorder - 1;
        end
    end
    K_i = length(hmm_i.state);

    if ii==1
        assig = 1:K_i;
        if K_i<K
            warning('The first HMM run needs to return K states, you might want to start again..\n')
        end
        hmm_init = struct('train',hmm_i.train);
        hmm_init.train.active = ones(1,K);
       for as=1:K_i
       hmm_init.state(as).Omega=hmm_i.state(as).Omega;
       hmm_init.state(as).W=hmm_i.state(as).W;
       end
    end

end
if isempty(options.BIGprior)
    for k = 1:K
        hmm_init.state(k).prior = hmm_i.state(k).prior;
        if isfield(hmm_init.state(k).prior,'Omega')
            if strcmp(options.covtype,'diag')
                hmm_init.state(k).prior.Omega.Gam_rate = 0.5 * range_data;
            elseif strcmp(options.covtype,'full')
                hmm_init.state(k).prior.Omega.Gam_rate = diag(range_data);
            end
        end
        if isfield(hmm_init.state(k).prior,'Mean')
            hmm_init.state(k).prior.Mean.Mu = zeros(ndim,1);
            hmm_init.state(k).prior.Mean.S = ((range_data/2).^2)';
            hmm_init.state(k).prior.Mean.iS = 1 ./ hmm_init.state(k).prior.Mean.S;
        end
    end
    hmm_init.prior = hmm_i.prior;
    if strcmp(options.covtype,'uniquediag')
        hmm_init.prior.Omega.Gam_rate = 0.5 * range_data;
    elseif strcmp(options.covtype,'uniquefull')
        hmm_init.prior.Omega.Gam_rate = diag(range_data);
    end
else
    for k = 1:K
        hmm_init.state(k).prior = options.BIGprior.state(k).prior;
    end
    hmm_init.prior.Dir2d_alpha = options.BIGprior.Dir2d_alpha;
    hmm_init.prior.Dir_alpha = options.BIGprior.Dir_alpha;
end
hmm_init.K = K;
hmm_init.train.BIGNbatch = options.BIGNbatch;
hmm_init.train.Sind = Sind; 

% update transition probabilities
hmm_init.Dir_alpha = sum(Dir_alpha_init,2)' + Dir_alpha_prior;
hmm_init.Dir2d_alpha = sum(Dir2d_alpha_init,3) + Dir2d_alpha_prior;
[hmm_init.P,hmm_init.Pi] =  computePandPi(hmm_init.Dir_alpha,hmm_init.Dir2d_alpha);

% Compute free energy
hmm_init_i = hmm_init;
for i = 1:N
    [X,XX,Y,Ti] = loadfile(Xin{i},T{i},options);
    data = struct('X',X,'C',NaN(size(XX,1),K));
    [Gamma,~,Xi,l] = hsinference(data,Ti,hmm_init_i,Y,[],XX);
    checkGamma(Gamma,Ti,hmm_init_i.train,i);
    subjfe_init(i,1:2) = evalfreeenergy([],Ti,Gamma,Xi,hmm_init_i,[],[],[1 0 1 0 0]); % Gamma entropy&LL
    loglik_init(i) = sum(l);
end
subjfe_init(:,3) = evalfreeenergy([],[],[],[],hmm_init,[],[],[0 0 0 1 0]) / N; % "share" P and Pi KL
statekl_init = sum(evalfreeenergy([],[],[],[],hmm_init,[],[],[0 0 0 0 1])); % state KL
fe = - sum(loglik_init) + sum(subjfe_init(:)) + statekl_init;

if fe<best_fe
    best_fe = fe;
    hmm = hmm_init;
    info.Dir2d_alpha = Dir2d_alpha_init; info.Dir_alpha = Dir_alpha_init;
    info.subjfe = subjfe_init;
    info.loglik = loglik_init;
    info.statekl = statekl_init;
    info.fehist = (-sum(info.loglik) + sum(info.statekl) + sum(sum(info.subjfe)));
end


hmm.prior.Dir_alpha_prior = Dir_alpha_prior;
hmm.prior.Dir2d_alpha_prior = Dir2d_alpha_prior;

