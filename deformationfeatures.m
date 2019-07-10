function [rdsv]=deformationfeatures(Ttilde,x,y,rankmin);
% Computation of the rate of decay of the SV curve
% INPUT
% Ttilde:  estimated transformation
% x:       patterns in the ROIX
% y:       patterns in the ROIY
% rankmin: minimum between the rank of x and the rank of y
% OUTPUT
% rdsv:    rate of decay of an exponential function fitted to the SV curve
% Alessio Basti 20/02/2019 (Basti et al. 2019)

[U K V]=svd(Ttilde);
f=fit((1:rankmin)',diag(K(1:rankmin,1:rankmin)),'exp1');
rdsv=f.b;

return