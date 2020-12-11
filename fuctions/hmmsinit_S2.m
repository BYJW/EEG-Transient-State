function [hmmi,Gamma,Xi] = hmmsinit_S2(data,T,options,h)

%%%%% to do a initialization for the model training
removesoptions;
options = rmfield(options,'orders');
options.pca = 0; % done in loadfile.m
options.pca_spatial = 0;
options.embeddedlags = 0; % done in loadfile.m
options.filter = [];
options.detrend = 0; 
options.onpower = 0; 
options.downsample = 0; % done in loadfile.m
options.grouping = [];
options.initTestSmallerK = 0;
options.BIGbase_weights = 0.95;
Sind = options.Sind;

N = length(T);
if iscell(T)
    for i = 1:length(T)
        if size(T{i},1)==1, T{i} = T{i}'; end
    end
    if size(T,1)==1, T = T'; end
    T = cell2mat(T);
end
checkdatacell;
%%%%%%% to do the initialization
[options,data] = checkoptions(options,data,T,0);
data = standardisedata(data,T,options.standardise); 
hmm_wr = struct('train',struct());
hmm_wr.K = options.K;

hmm_wr.train = options;
hmm_wr = hmmhsinit(hmm_wr);  %%%%%%% ****** intiate the hmm.Dir2d_alpha_prior  hmm.P
%%%%%%%%  measure the residuals
[residuals0,W0] =  getresiduals(data.X,T,hmm_wr.train.Sind,hmm_wr.train.maxorder,hmm_wr.train.order,...
    hmm_wr.train.orderoffset,hmm_wr.train.timelag,hmm_wr.train.exptimelag,hmm_wr.train.zeromean);
%%%%%%%%  use the current loaded hmm to be model for training
hmm_wr.state=h.hmm.state;  

%%%%% the hsinference function to get the Gamma and update hmm.Pi hmm.P
hmm0 = hmm_wr;
hmm0.train.active = ones(1,hmm0.K);
[hmmi,Gamma,Xi,fehist] = hmmtrain_S(data,T,hmm0,[],residuals0,[]); 
