function [rdsv]=deformationfeatures(Ttilde,x,y);
% Computation of the rate of decay of the SV curve
% INPUT
% Ttilde: estimated transformation
% x:      patterns in the ROIX
% y:      patterns in the ROIY
% OUTPUT
% rdsv:   rate of decay of an exponential function fitted to the SV curve
% Alessio Basti 20/02/2019 (Basti et al. 2019)

[U K V]=svd(Ttilde);
f=fit((1:size(x,2))',diag(K(1:size(x,2),1:size(x,2))),'exp1');
rdsv=f.b;

return