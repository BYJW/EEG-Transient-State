function [ res ] = gauss_model_2fit( data , plotflag )



res=[];

data=squash(data);

normalised_logistic_data=normalise(data);   
xx=min(normalised_logistic_data):range(normalised_logistic_data)/50:max(normalised_logistic_data);


   
mix = gmm(1, 2, 'diag');
mix = gmminit(mix, normalised_logistic_data, foptions);
mix = gmmem(mix, normalised_logistic_data, foptions);
[prob] = gmmactiv(mix, xx');
prob2=prob;
for cc=1:length(mix.centres)
    prob2(:,cc)=prob2(:,cc)*mix.priors(cc); 
end

[mix.centres,ort]=sort(mix.centres,'ascend');
prob2=prob2(:,ort);
cen(1)=nearest(prob2(:,end-1),max(prob2(:,end-1)));
cen(2)=nearest(prob2(:,end),max(prob2(:,end)));

th=xx(nearest(prob2(cen(1):cen(2),end-1)-prob2(cen(1):cen(2),end),0)+cen(1)-1);
orig_th=data(nearest(normalised_logistic_data,th));



%%%%%%%%%%%%%%%%%%%%%
%% do plots

%%
if plotflag
figure;
set(gcf,'Position',[440   378   1000   420])

subplot(121);
plot(xx,prob2,'LineWidth',3)
[prob3] = gmmprob(mix, xx');
plot4paper('normalised edge strength','');

subplot(122);
hh=hist(normalised_logistic_data,xx);
bar(xx,hh/sum(hh));
ho;
plot(xx,prob3/sum(prob3),'r','LineWidth',3);
title(num2str(th));
plot4paper('normalised edge strength','');
    
end

res.mix=mix;
res.normalised_th=th;
res.orig_th=orig_th;
res.data=normalised_logistic_data;

end

