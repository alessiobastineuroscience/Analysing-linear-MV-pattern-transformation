function [T]=simulmatrixspecificSV(dimx,dimy,deformation); 
% Starting from a matrix whose entries are realizations of Gaussian random
% variable, this function changes the singular values in order to allow them 
% to have a specific rate of decay.
% INPUT
% dimx:          number of voxels in ROIX
% dimy:          number of voxels in ROIY
% deformation:   rate of decay
% OUTPUT
% T:             simulated transformation with a specific decay of the singular values
% Alessio Basti 20/02/2019 (Basti et al. 2019)

T=randn(dimy,dimx);
[U K V]=svd(T);
eigsim=exp(-deformation*(1:min(size(T))));
K(1:(min(size(T))),1:min(size(T)))=diag(eigsim);
T=U*K*V';

return