function [hmm,Gamma,Xi,fehist] = hmmtrain_S(data,T,hmm,Gamma,residuals,fehist)
%
% Train Hidden Markov Model using using Variational Framework


if nargin<6, fehist=[]; end
cyc_to_go = 0;
setxx;

for cycle=1:hmm.train.cyc
    
    if hmm.train.updateGamma
        
        %%%% E step
        if hmm.K>1 || cycle==1
            % state inference
            [Gamma,~,Xi] = hsinference(data,T,hmm,residuals,[],XX);
            status = checkGamma(Gamma,T,hmm.train);
            % check local minima
            epsilon = 1;
            while status == 1
                hmm = hmmperturb(hmm,epsilon); %%%%%%******** generate the hmm.Omega.Gam_rate, shape
                disp('Stuck in bad local minima - perturbing the model and retrying...')
                [Gamma,~,Xi] = hsinference(data,T,hmm,residuals,[],XX);
                status = checkGamma(Gamma,T,hmm.train);
                epsilon = epsilon * 2;
            end
%             
%             % any state to remove?
            [as,hmm,Gamma,Xi] = getactivestates(hmm,Gamma,Xi);
            if hmm.train.dropstates
                if any(as==0)
                    cyc_to_go = hmm.train.cycstogoafterevent;
                    data.C = data.C(:,as==1);
                    [Gamma,~,Xi] = hsinference(data,T,hmm,residuals,[],XX);
                    checkGamma(Gamma,T,hmm.train);
                end
                if sum(hmm.train.active)==1
                    fehist(end+1) = sum(evalfreeenergy(data.X,T,Gamma,Xi,hmm,residuals,XX));
                    if hmm.train.verbose
                        fprintf('cycle %i: All the points collapsed in one state, free energy = %g \n',...
                            cycle,fehist(end));
                    end
                    K = 1; break
                end
            end
            setxx;
        end

        %%%% Free energy computation
        fehist(end+1) = sum(evalfreeenergy(data.X,T,Gamma,Xi,hmm,residuals,XX));
        strwin = ''; if hmm.train.meancycstop>1, strwin = 'windowed'; end
        if cycle>(hmm.train.meancycstop+1) 
            chgFrEn = mean( fehist(end:-1:(end-hmm.train.meancycstop+1)) - ...
                fehist(end-1:-1:(end-hmm.train.meancycstop)) )  ...
                / (fehist(1) - fehist(end));
            if hmm.train.verbose
                fprintf('cycle %i free energy = %g, %s relative change = %g \n',...
                    cycle,fehist(end),strwin,chgFrEn); 
            end
            if (abs(chgFrEn) < hmm.train.tol) && cyc_to_go==0, break; end
        elseif hmm.train.verbose
            fprintf('cycle %i free energy = %g \n',cycle,fehist(end)); %&& cycle>1
        end
        if cyc_to_go>0, cyc_to_go = cyc_to_go - 1; end
        
    else
        Xi=[]; fehist=0;
    end
     
        hmm = hsupdate(Xi,Gamma,T,hmm);    %%%%%%******** finally updatae the transition and initial state probability

end

for k = 1:K
    if isfield(hmm.state(k),'cache')
        hmm.state(k) = rmfield(hmm.state(k),'cache');
    end
end

if hmm.train.verbose
    fprintf('Model: %d states, %d data samples, covariance: %s \n', ...
        K,sum(T),hmm.train.covtype);
    if hmm.train.exptimelag>1
        fprintf('Exponential lapse: %g, order %g, offset %g \n', ...
            hmm.train.exptimelag,hmm.train.order,hmm.train.orderoffset)
    else
        fprintf('Lapse: %d, order %g, offset %g \n', ...
            hmm.train.timelag,hmm.train.order,hmm.train.orderoffset)
    end
    if hmm.train.useMEX==0
        fprintf('MEX file was not used \n')
    else
        fprintf('MEX file was used for acceleration \n')
    end
end

end
