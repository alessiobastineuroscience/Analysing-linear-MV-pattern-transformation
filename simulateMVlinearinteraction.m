function [x,y]=simulateMVlinearinteraction(dimx,dimy,dimt,levelofnoise,levelofsparsity,levelofdeformation,numberofsubjs);
% Simulate linear MV-interaction between ROIX and ROIY
% INPUT
% dimx:               number of voxels in the first ROI (ROIX)
% dimy:               number of voxels in the second ROI (ROIY)
% dimt:               number of stimuli
% levelofnoise:       weight of the noise in the simulated interaction
% levelofsparsity:    percentage of sparsity of the transformation
% levelofdeformation: simulated rate of decay of singular value curve.
%                     Either levelofsparsity or levelofdeformation has to be set to [].
% numberofsubjs:      number of subjects
% OUTPUT
% x:                  each cell represents the matrix with demeaned and standardized MV-pattern   
%                     for the ROIX for a specific subject and run
% y:                  each cell represtns the MV-patterns for the ROIY for a specific
%                     subject and run
% Alessio Basti 20/02/2019 (Basti et al. 2019)

if(length(levelofsparsity)==1)
    for nsubj=1:numberofsubjs
        for mruns=1:2
            x{mruns,nsubj}=randn(dimx,dimt);
            for i=1:size(x{mruns,nsubj},2)
                x{mruns,nsubj}(:,i)=(x{mruns,nsubj}(:,i)-mean(x{mruns,nsubj}(:,i)))/std(x{mruns,nsubj}(:,i));
            end
            % simulate a transformation with a specific percentage of sparsity
            T=simulsparsematrix(dimx,dimy,levelofsparsity);
            noise=randn(dimy,dimt);
            y{mruns,nsubj}=(1-levelofnoise)*T*x{mruns,nsubj}/norm(T*x{mruns,nsubj},'fro')+levelofnoise*noise/norm(noise,'fro');
            for i=1:size(y{mruns,nsubj},2)
                y{mruns,nsubj}(:,i)=(y{mruns,nsubj}(:,i)-mean(y{mruns,nsubj}(:,i)))/std(y{mruns,nsubj}(:,i));
            end
        end
    end
else
    for nsubj=1:numberofsubjs
        for mruns=1:2
            x{mruns,nsubj}=randn(dimx,dimt);
            for i=1:size(x{mruns,nsubj},2)
                x{mruns,nsubj}(:,i)=(x{mruns,nsubj}(:,i)-mean(x{mruns,nsubj}(:,i)))/std(x{mruns,nsubj}(:,i));
            end
            % simulate a transformation with a specific decay of SVs
            T=simulmatrixspecificSV(dimx,dimy,levelofdeformation);  
            noise=randn(dimy,dimt);
            y{mruns,nsubj}=(1-levelofnoise)*T*x{mruns,nsubj}/norm(T*x{mruns,nsubj},'fro')+levelofnoise*noise/norm(noise,'fro');
            for i=1:size(y{mruns,nsubj},2)
                y{mruns,nsubj}(:,i)=(y{mruns,nsubj}(:,i)-mean(y{mruns,nsubj}(:,i)))/std(y{mruns,nsubj}(:,i));
            end
        end
    end

end