clear;clc;close all
% Starting from two matrices x and y composed of demeaned and standardized MV-patterns 
% (either simulated or real), this script can be used for estimating 
% 1) the linear transformation between x and y by using the ridge regression method, 
% 2) the goodness-of-fit (GOF) and 3) the induced pattern deformation 
% (via Monte Carlo procedure which takes into account both the GOF value and 
% the rate of decay of the singular value curve(RDSV)) as in Basti et al. 2019.

%%   
% In order to understand how the script works, let us simulate the data 
dimx=100; % number of voxels in the first ROI (ROIX)
dimy=150; % number of voxels in the second ROI (ROIY)
dimt=90;  % number of stimuli
levelofnoise=0.3; % weight of the noise, the one of the signal is (1-levelofnoise)
levelofsparsity=[]; % in this toy example let us only consider the deformation
levelofdeformation=0.5; % 0 denotes an orthogonal transformation, the higher
                        % the parameter, the higher the deformation induced 
                        % by the transformation
numberofsubjs=4; % number of subjects to simulate
% simulate linear MV-interaction between ROIX and ROIY
[x,y]=simulateMVlinearinteraction(dimx,dimy,dimt,levelofnoise,levelofsparsity,levelofdeformation,numberofsubjs);

%%
% Let us estimate the transformation and the metrics
lambdas=10.^(-2:0.1:5); %set of regularization parameters
results=featuresevaluation(x,y,lambdas);

% In order to associate to each pair (GOF, RDSV) a specific level of
% deformation, let us use the Monte Carlo (MC) approach
deformations=[0, 0.01, 0.1, 1]; % rate of decay of singular values to simulate
numberofrep=10; % number of simulations for each level of deformation and noise in the MC
approach=0; %using lambdas in the set. Instead approach=1 means using the lambdas obtained before
reglambda.optimal=results.lambdas;
reglambda.set=lambdas;
resultsMCdeformation=MCproceduredeformation(x,dimx,dimy,dimt,reglambda,approach,deformations,results.rank,numberofrep);

%% 
% plot (equivalent to that used in Basti et al. 2019). The black square (
% with the error bars) denotes the average (and std) GOF and RDSV across the
% subjects (either on real or on the data simulated as above), while the 
% coloured squares denotes the results obtained with the MC approach. By
% looking at the position of the black square with respect to those of the
% coloured squares, it is possible to estimate the level of deformation
% induced by the transformation of interest: for instance, if the
% black square lies between the curve representing the orthogonal transformation
% (i.e. 0 deformation) and  of 
% sparsity in the MC, it means that the estimated level of deformation
% for the transformation of interest is in the range 0-0.01.
color={'m';'c';'g';'r';'b';'y';'k'};
figure('Position',[0 0 900 900])
normfactor=sqrt(numberofrep*numberofsubjs*2);
normfactorsubj=sqrt(numberofsubjs*2);
for idef=1:numel(deformations)
   hold on
   errorbar(mean(squeeze(resultsMCdeformation.gof(idef,:,:)),2),mean(squeeze(resultsMCdeformation.rdsv(idef,:,:)),2),std(squeeze(resultsMCdeformation.rdsv(idef,:,:)),0,2)/normfactor,std(squeeze(resultsMCdeformation.rdsv(idef,:,:)),0,2)/normfactor,std(squeeze(resultsMCdeformation.gof(idef,:,:)),0,2)/normfactor,std(squeeze(resultsMCdeformation.gof(idef,:,:)),0,2)/normfactor,'-sb','MarkerSize',6,'MarkerEdgeColor',color{idef},'MarkerFaceColor',color{idef},'DisplayName',strcat('Rate of decay SV=',num2str(deformations(idef))))
end
hold on
set(gca,'FontSize',20)
errorbar(mean(results.gof),mean(results.rdsv),std(results.rdsv)/normfactorsubj,std(results.rdsv)/normfactorsubj,std(results.gof)/normfactorsubj,std(results.gof)/normfactorsubj,'-sk','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor','k','DisplayName','Estimate on data of interest')
title(strcat('Estimated level of deformation (original value=',num2str(levelofdeformation),')'))
xlabel('Goodness-of-fit')
ylabel('Rate of decay of the singular values')
legend('Location','southeast')
axis([0 100 2*mean(results.rdsv) 0])

%save('Results_deformation.mat','results','resultsMCdeformation') 