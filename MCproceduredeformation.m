function [resultsMCdeformation]=MCproceduredeformation(x,dimx,dimy,dimt,reglambda,approach,deformations,rank,numberofrep);
% Monte Carlo procedure for pattern deformation.
% INPUT
% x:                    MV-pattern of interest for the input region
% dimx:                 number of voxels in the ROIX
% dimt:                 number of stimuli
% reglambda:            regularization parameter
% approach:             1 means using the optimal lambda obtained on real
%                       data, while 0 means using lambdas in the set given as input
% deformations:         different rates of decay of SVs to investigate
% rank:                 minimum rank of the real data
% numberofrep:          number of simulations for each investigated level of
%                       deformation and noise
% OUTPUT
% resultsMCdeformation: structure with goodness-of-fit, the rate of decay
%                       of SVs
% Alessio Basti 20/02/2019 (Basti et al. 2019)

gammas=0:0.1:0.9;
totalns=numberofrep*numel(gammas)*numel(deformations)*size(x,2)*2;
count=0;

for mruns=1:2
    for nsubj=1:size(x,2)
        X=x{mruns,nsubj};
        for jrep=1:numberofrep
            for kgam=1:numel(gammas)
                for idef=1:numel(deformations)
                    clearvars -except x X gammas deformations reglambda approach indexesr indexesc dimt dimx dimy numberofrep idef kgam jrep mruns nsubj resultsMCdeformation rank count totalns

                    % simulate transformation and MV-patterns of the ROIY
                    T=simulmatrixspecificSV(dimx,dimy,deformations(idef));         
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

                    % computation of the rate of decay of the SV curve
                    [rdsv]=deformationfeatures(Ttilde,X,y,rank(nsubj,mruns));

                    % create a structure cointaining: the goodness-of-fit, the
                    % density curve and its rate of decay
                    resultsMCdeformation.gof(idef,kgam,jrep+numberofrep*(nsubj-1)+(mruns-1)*numberofrep*size(x,2))=gof;
                    resultsMCdeformation.rdsv(idef,kgam,jrep+numberofrep*(nsubj-1)+(mruns-1)*numberofrep*size(x,2))=rdsv;
                   
                    count=count+1;
                    100*(count)/totalns
                end
            end
        end
    end
end
resultsMCdeformation.analysedlevofdeformation=deformations;
return