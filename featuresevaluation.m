function [results]=featuresevaluation(x,y,lambdas);
% Estimate of the transformations and the metrics for data of interest.
% INPUT
% x:         fMRI MV-pattern for the input region (each cell denotes a 
%            run for a specific subj)
% y:         fMRI MV-pattern for the output region (each cell denotes a
%            run for a specific subj)
% lambdas:   set of possible regulariation parameter
% OUTPUT
% results:   structure with estimated transformation, goodness-of-fit (GOF), 
%            density curve its rate of decay
% Alessio Basti 20/02/2019 (Basti et al. 2019)

% in case of real data it is better to use an across-sessions approach  
for nsubj=1:size(x,2)
    for mruns=1:2       
        X=x{mruns,nsubj};
        Y=y{mruns,nsubj};
        
        % computation of optimal parameter, of the estimated transformation and GOF
        [estT,optlambda(nsubj,mruns),gof(mruns)]=ridgeregmethod(X,Y,lambdas);
        
        % computation of the density curve and of its rate of decay
        [rdd(mruns),density(:,mruns)]=sparsityfeatures(estT,X,Y);
        
        % computation of the rate of decay of the SV curve
        rankmin(nsubj,mruns)=min([rank(X),rank(Y)]);
        [rdsv(mruns)]=deformationfeatures(estT,X,Y,rankmin(nsubj,mruns));
        
        Ttilde(mruns,1:size(Y,1),1:size(X,1))=estT;
    end       
    % structure cointaining the results
    results.estimtransf{nsubj}=Ttilde;
    results.gof(nsubj)=mean(gof);
    results.rdd(nsubj)=mean(rdd);
    results.density(:,nsubj)=mean(density,2);
    results.rdsv(nsubj)=mean(rdsv);    
end
results.lambdas=optlambda;
results.rank=rankmin;

return