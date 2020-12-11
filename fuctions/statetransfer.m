function [Gamma, vpath, hmm]=statetransfer(data,T,options,h)


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

