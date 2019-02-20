function [T]=simulsparsematrix(dimx,dimy,sparsity);
% Starting from a matrix whose entries are realizations of Gaussian random
% variable, this function removes a percentage of components in order to allow 
% the matrix to have a specific percentage of sparsity.
% INPUT
% dimx:       number of voxels in ROIX
% dimy:       number of voxels in ROIY
% sparsity:   percentage of sparsity
% OUTPUT
% T:          simulated transformation
% Alessio Basti 20/02/2019 (Basti et al. 2019)

indexesr=randi(dimy,1,250000);
indexesc=randi(dimx,1,250000);

degreeofsparsity=sparsity/100;
T=randn(dimy,dimx);
N=floor(prod(size(T))*degreeofsparsity);
i=1;
k=1;
while i<N+1
    if(abs(T(indexesr(k),indexesc(k)))>0)
        T(indexesr(k),indexesc(k))=0;
        i=i+1;
    end
    k=k+1;
end

return