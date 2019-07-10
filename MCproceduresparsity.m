function [resultsMCsparsity]=MCproceduresparsity(x,dimx,dimy,dimt,reglambda,approach,sparsities,numberofrep);
% Monte Carlo procedure for sparsity.
% INPUT
% x:                 MV-pattern matrix of interest for the input region
% dimx:              number of voxels in the ROIX
% dimt:              number of stimuli
% reglambda:         regulariation parameters
% approach:          1 means using the optimal lambda obtained on real
%                    data, while 0 means using lambdas in the set given as input
% sparsities:        percentages of sparsity to simulate
% numberofrep:       number of simulations for each investigated level of
%                    sparsity and noise
% OUTPUT
% resultsMCsparsity: structure with goodness-of-fit, density curve and its rate of decay
% Alessio Basti 20/02/2019 (Basti et al. 2019)

gammas=0:0.1:0.9;
totalns=numberofrep*numel(gammas)*numel(sparsities)*size(x,2)*2;
count=0;

for mruns=1:2
    for nsubj=1:size(x,2)
        X=x{mruns,nsubj};
        for jrep=1:numberofrep
            for kgam=1:numel(gammas)
                for ispar=1:numel(sparsities)
                    clearvars -except x X gammas sparsities reglambda approach indexesr indexesc dimt dimx dimy numberofrep ispar kgam jrep mruns nsubj resultsMCsparsity count totalns

                    % simulate transformation and MV-patterns of the ROIY
                    T=simulsparsematrix(dimx,dimy,sparsities(ispar));         
                    noise=randn(dimy,dimt);
                    y=(1-gammas(kgam))*T*X/norm(T*X,'fro')+gammas(kgam)*noise/norm(noise,'fro');
                    for i=1:size(y,2)
                        y(:,i)=(y(:,i)-mean(y(:,i)))/std(y(:,i));
                    end

                    % computation of the estimated transformation and GOF
                    if(approach==1)
                        [Ttilde,optlambda,gof]=ridgeregmethod(X,y,reglambda.optimal(nsubj,mruns));
                    elseif(approach==0)
                        [Ttilde,optlambda,gof]=ridgeregmethod(X,y,reglambda.set);
                    end
                    % computation of the density curve and of its rate of
                    % decay
                    [rdd,density]=sparsityfeatures(Ttilde,X,y);

                    % create a structure cointaining: the goodness-of-fit, the
                    % density curve and its rate of decay
                    resultsMCsparsity.gof(ispar,kgam,1:100,jrep+numberofrep*(nsubj-1)+(mruns-1)*numberofrep*size(x,2))=gof;
                    resultsMCsparsity.rdd(ispar,kgam,jrep+numberofrep*(nsubj-1)+(mruns-1)*numberofrep*size(x,2))=rdd;
                    resultsMCsparsity.density(ispar,kgam,1:100,jrep+numberofrep*(nsubj-1)+(mruns-1)*numberofrep*size(x,2))=density;

                    count=count+1;
                    100*(count)/totalns
                end
            end
        end
    end
end
resultsMCsparsity.analysedpercofsparsity=sparsities;
return