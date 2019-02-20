clear;clc;close all

%% either the real or the simulated fMRI MV-patterns for the input ROIX.
dimx=100; %number of voxels in the first ROI (ROIX)
dimy=150; %number of voxels in the second ROI (ROIY)
dimt=90;  %number of stimuli
levelofnoise=0.4;
levelofdeformation=[];
levelofsparsity=75; %percentage of sparsity of the transformation
numberofsubjs=4; %number of subjects
% simulate linear MV-interaction between ROIX and ROIY
[x,y]=simulateMVlinearinteraction(dimx,dimy,dimt,levelofnoise,levelofsparsity,levelofdeformation,numberofsubjs);

%%
lambdas=10.^(-2:0.1:5); %set of possible regularization parameter
sparsities=[50,70,80,90]; % percentages of sparsity to simulate
gammas=0:0.1:0.9;  % the strength of noise in the model
numberofrep=1; %number of simulations for each level of sparsity and gamma

results=featuresevaluation(x,y,lambdas);
resultsMCsparsity=MCproceduresparsity(x,dimx,dimy,dimt,lambdas,sparsities,gammas,numberofrep);

%save('Results_sparsity.mat','results','resultsMCsparsity','sparsities') 

%% plot (equivalent to that ..)
color={'m';'c';'g';'r';'b'};
figure('Position',[50 50 900 900])
for ispar=1:numel(sparsities)
   hold on
   errorbar(mean(squeeze(resultsMCsparsity.gof(ispar,:,1,:)),2),mean(squeeze(resultsMCsparsity.rdd(ispar,:,:)),2),std(squeeze(resultsMCsparsity.rdd(ispar,:,:)),0,2),std(squeeze(resultsMCsparsity.rdd(ispar,:,:)),0,2),std(squeeze(resultsMCsparsity.gof(ispar,:,1,:)),0,2),std(squeeze(resultsMCsparsity.gof(ispar,:,1,:)),0,2),'-sb','MarkerSize',6,'MarkerEdgeColor',color{ispar},'MarkerFaceColor',color{ispar},'DisplayName',strcat('Sparsity=',num2str(sparsities(ispar)),'%'))
   %legend(strcat('Sparsity=',num2str(sparsities(ispar))))
end
hold on
errorbar(mean(results.gof),mean(results.rdd),std(results.rdd),std(results.rdd),std(results.gof),std(results.gof),'-sk','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName','Estimate on real data')
title('Estimate of sparsity')
xlabel('Goodness-of-fit')
ylabel('Rate of decay of the density')
legend('Location','southwest')

