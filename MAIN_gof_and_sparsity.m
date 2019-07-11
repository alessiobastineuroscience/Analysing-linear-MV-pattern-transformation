clear;clc;close all
% Starting from two matrices x and y composed of demeaned and standardized MV-patterns 
% (either simulated or real), this script can be used for estimating 
% 1) the linear transformation between x and y by using the ridge regression method, 
% 2) the goodness-of-fit (GOF) metric and 3) the percentage of sparsity
% (via Monte Carlo procedure that takes into account both the GOF
% value and the rate of decay of the density curve (RDD)) as in Basti et al. 2019.

%%   
% In order to understand how the script works, let us simulate the data 
dimx=100; % number of voxels in the first ROI (ROIX)
dimy=150; % number of voxels in the second ROI (ROIY)
dimt=90;  % number of stimuli
levelofnoise=0.3; % weight of the noise, the one of the signal is thus (1-levelofnoise)
levelofsparsity=65; % simulated percentage of sparsity. 0 denotes two fully 
                    % connected regions while e.g. 80 denotes a transformation
                    % in which the 80% of entries are equal to 0, i.e. on average
                    % each voxel in the ROIX interact with only 20%
                    % of the voxels in the other region
levelofdeformation=[]; % let us only consider the sparsity
numberofsubjs=4; % number of subjects to simulate
% simulate linear MV-interaction between ROIX and ROIY
[x,y]=simulateMVlinearinteraction(dimx,dimy,dimt,levelofnoise,levelofsparsity,levelofdeformation,numberofsubjs);

%%
% Let us estimate now the transformations and the metrics
lambdas=10.^(-2:0.1:5); %set of regularization parameters
results=featuresevaluation(x,y,lambdas);

% In order to associate to each pair (GOF, RDD) a specific percentage of 
% sparsity, let us use the Monte Carlo (MC) approach
sparsities=[0,50,60,70,80]; % percentage of sparsity to investigate
numberofrep=10; % number of repetitions for each level of sparsity and noise in the MC
approach=0; %using lambdas in the set. Instead approach=1 means using the lambdas obtained before
reglambda.optimal=results.lambdas;
reglambda.set=lambdas;
resultsMCsparsity=MCproceduresparsity(x,dimx,dimy,dimt,reglambda,approach,sparsities,numberofrep);

%% 
% plot (equivalent to that used in Basti et al. 2019). The black square
% (with the error bars) denotes the mean (and standard error of the mean) GOF 
% and RDD across the subjects (either on real or on the data simulated as above),  
% while the coloured squares denotes the results obtained with the MC approach. 
% By looking at the position of the black square with respect to those of the
% coloured squares, it is possible to estimate the percentage of sparsity
% of the transformation of interest: for instance, if the
% black square lies between the curve representing the 70% and the 80% of 
% sparsity in the MC, it means that the estimated percentage of sparsity 
% for the transformation of interest is in the range 70-80%.
color={'m';'c';'g';'r';'b';'y';'k'};
figure('Position',[0 0 900 900])
normfactor=sqrt(numberofrep*numberofsubjs*2);
normfactorsubj=sqrt(numberofsubjs*2);
for ispar=1:numel(sparsities)
   hold on
   errorbar(mean(squeeze(resultsMCsparsity.gof(ispar,:,1,:)),2),mean(squeeze(resultsMCsparsity.rdd(ispar,:,:)),2),std(squeeze(resultsMCsparsity.rdd(ispar,:,:)),0,2)/normfactor,std(squeeze(resultsMCsparsity.rdd(ispar,:,:)),0,2)/normfactor,std(squeeze(resultsMCsparsity.gof(ispar,:,1,:)),0,2)/normfactor,std(squeeze(resultsMCsparsity.gof(ispar,:,1,:)),0,2)/normfactor,'-sb','MarkerSize',10,'MarkerEdgeColor',color{ispar},'MarkerFaceColor',color{ispar},'DisplayName',strcat('Sparsity=',num2str(sparsities(ispar)),'%'),'LineWidth',1.5)
end
hold on
set(gca,'FontSize',20)
errorbar(mean(results.gof),mean(results.rdd),std(results.rdd)/normfactorsubj,std(results.rdd)/normfactorsubj,std(results.gof)/normfactorsubj,std(results.gof)/normfactorsubj,'-sk','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName','Estimate for the data of interest')
title(strcat('Estimated percentage of sparsity (original value=',num2str(levelofsparsity),'%)'))
xlabel('Goodness-of-fit')
ylabel('Rate of decay of the density')
legend('Location','southwest')
axis([0 100 2*mean(results.rdd) 0])

%save('Results_sparsity.mat','results','resultsMCsparsity') 
