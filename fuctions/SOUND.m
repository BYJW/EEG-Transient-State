function [corrected_data, sigmas, lambda, x, i_ref, k] = SOUND(data, LFM,rf, lambda, i_ref, max_iter)
%
%This function performs SOUND algorithm for a given data
%
%Input
% data: 2D or 3D array of channels x times x trials(optional)
% LFM: lead-field matrix
% rf: small number e.g. 1e-4 to provide a regularization factor for
% computing the initial data-driven estimation. Increase this if Matlab
% warns about ill-conditioned matrix.
% lambda: regularization parameter ~ 1/SNR. If lambda=[], it'll be
% automatically adjusted
% i_ref: reference channel, with as little noise as possible. 
% If iref = 0, the current reference is kept 
% (Note that LFM must then have the same reference!). If iref=[], it'll be
% automatically adjusted.
% max_iter: maximum number of iterations
%
%Output
% corrected_data: in 2D channels x times in the given/automatically selected reference
% sigmas: found noise variaces over channels
% lambda: the final regualization parameter
% x: minimum-norm-estimated sources
% i_ref: the used reference channel
% k: iterations untill convergence

[M,T,rep]=size(data);

if rep> 1
    data2=reshape(data, M, []);
    [y_solved, sigmas] = simple_wiener(mean(data2,3), rf);
    y_solved=reshape(y_solved, [M, T, rep]);
   
    data=mean(data,3);
    y_solved=mean(y_solved,3);
    sigmas=sigmas./rep;
else 
    [y_solved, sigmas] = simple_wiener(data,rf);
end

if isempty(i_ref)
    [~ , i_ref]=min(sigmas);
end

if i_ref ==0    
    data_rest=data;
    LFM_rest=LFM;
    y_solved_rest=y_solved;    
    noise_cov = diag(sigmas); 
    Mrest=M;
else
    rest_ch=setdiff(1:M,i_ref);
    data_rest=data(rest_ch,:)-repmat(data(i_ref,:),M-1,1);
    LFM_rest=LFM(rest_ch,:)-repmat(LFM(i_ref,:), M-1,1);
    y_solved_rest=y_solved(rest_ch,:)-repmat(y_solved(i_ref,:), M-1, 1);
        
    noise_cov = diag(sigmas+sigmas(i_ref)); %!!!
    noise_cov=(noise_cov(rest_ch, rest_ch));
    
    Mrest=M-1;
end
sigmas_old=diag(noise_cov);
if isempty(lambda)
    lambda = trace(LFM_rest*(LFM_rest)')/trace(y_solved_rest*(y_solved_rest)'/T);
end
eps=1;
k=0;
while eps>1e-3 & k<max_iter
    k
    for i=1:Mrest
        chan = setdiff(1:Mrest,i);
        LFM_tmp = LFM_rest(chan,:);
        x = (LFM_tmp)'*((LFM_tmp*(LFM_tmp)' + lambda*noise_cov(chan,chan))\data_rest(chan,:));     
        y_solved_new(i,:) = LFM_rest(i,:)*x;
    end
    n = (y_solved_new-data_rest);
    sigmas_new=diag(n*n'./(T+1));
    noise_cov=diag(sigmas_new);%+noise_cov-diag(diag(noise_cov));
    y_solved_rest=y_solved_new;
    
    if isempty(lambda)
        lambda = trace(LFM_rest*(LFM_rest)')/trace(y_solved_rest*(y_solved_rest)'/T);
    end
    
    eps=max(abs((sigmas_new-sigmas_old)./sigmas_old))
    sigmas_old=sigmas_new;
    k=k+1;
end

sigmas=diag(noise_cov);
%k
x = (LFM_rest')*((LFM_rest*LFM_rest' + lambda*noise_cov)\data_rest);     
corrected_data = (LFM_rest)*x;


function [y_solved, sigmas] = simple_wiener(data,rf)
    [Nc]=size(data,1);
    C = cov(data');
    for i=1:Nc
        i
        idiff = setdiff(1:Nc,i);
        y_solved(i,:) = C(i,idiff)*((C(idiff,idiff)+rf/(Nc-1)*eye(Nc-1))\data(idiff,:));
    end
    sigmas = (diag((data-y_solved)*(data-y_solved)'))/(size(data,2));
end


end





