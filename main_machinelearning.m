clear all
clc
close all


%%  use the statistic test for feature selection
features=[features_temporal,features_psd,features_connectivity];
for i=1:size(features,2)
[h,p(i)]=ttest2(features(1:numVS,i),features(numVS+1:size(features,1),i));
end
refeatures=features(:,find(p<0.01));
feat=find(p<0.01);
 
%%   combined nested cross-validation with CARS-PLS-LDA model to classify the patients based on VS and MCS labels
addpath('/PLStoolbox');
precentest=10;  %%%% 10% for test, 90% for training
groups = mod(1:size(refeatures,1),10);
X=refeatures;
y=ones(1,size(refeatures,1));
y([length(y)-length(VSCRS)+1:end])=-1*ones(1,length(VSCRS));

loopnum=zeros(1,size(features,2));
scoreresult=zeros(1,size(refeatures,1));
accresult=zeros(1,precentest);
for group=0:precentest-1
    
testk = find(groups==group);  
calk = find(groups~=group);
Xcal=X(calk,:);ycal=y(calk);
Xtest=X(testk,:);ytest=y(testk);
    
A=6;
K=size(Xcal,2);  %%% leave one out
% K=3; %%% 4 fold cross-validation
method='center';
%%%% variables selection and parameter optimization
CARS1=carsplslda(Xcal,ycal',A,K,50,method,0); % simplified version with random elements removed so that results can be exactly reproducible
%%%%% the finaly selected variabiles
loopnum(feat(CARS1.vsel))=loopnum(feat(CARS1.vsel))+ones(1,length(feat(CARS1.vsel)));
%%%%% build the PLS-LDA model using the selected variables and the parameter
XX=Xcal(:,CARS1.vsel);
nLV=CARS1.optLV;
LDA=plslda(XX,ycal',nLV);
%%%%%  test the model using the testdata
Xtest=Xtest(:,CARS1.vsel);
Xtest=pretreat(Xtest,method,LDA.scale_para(1,:),LDA.scale_para(2,:));
Ttest=Xtest*LDA.Wstar;
testresult=Ttest*LDA.regcoef_lda_lv(1:end-1)+LDA.regcoef_lda_lv(end);
scoreresult(testk)=testresult;
end
accresult=sum(sign(scoreresult)==y)/length(y); 
ROC=roccurve(scoreresult,y',0); %+++ ROC area.
AUC=ROC.AUC;
k1=find(y==1);
k2=find(y==-1);
y_est=sign(scoreresult);
sensi=1-sum(y_est(k1)~=sign(y(k1)))/length(k1);
speci=1-sum(y_est(k2)~=sign(y(k2)))/length(k2);
figure
bar(loopnum);


X=score_crs';
Y=scoreresult;
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
scatter(X',Y,'blacko','LineWidth',1.5);
hold on
pp=polyfit(X',Y,1);
hold on
yy=polyval(pp,X');
plot(X',yy,'LineWidth',2,'Color',[0 0 1]);
[r2 rmse] = rsquare(Y,yy); %%%% used for regresssion
[r,p]=corr(score_crs',Y','type','Pearson');
r
p
xlabel('CRS-R');
xlim([0 18]);
ylim([-2.5 2.5]);
ylabel('Predicted scores');
set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2,'XTick',[0 5 10 15 20],'YTick',[-2 -1 0 1 2]); 


figure1 = figure;
axes1 = axes('Parent',figure1);
plot(scoreresult,'o','LineWidth',2);hold on
plot([1,65],[0 0],'color',[1 0 1],'LineWidth',2);
plot([25.5 25.5],[-2.5 2.5],'color',[0 1 1],'LineWidth',2);
% xlabel('CRS-R');
xlim([0 65]);
ylim([-2.5 2.5]);
ylabel('Predicted scores');
set(axes1,'FontSize',20,'FontWeight','bold','LineWidth',2,'XTick',[],'YTick',[-2 -1 0 1 2]); 

figure
[means,diffs,meanDiff,CR,linFit] = BlandAltman(score_crs, Y, 3);
figure
roccurve(scoreresult,y');



