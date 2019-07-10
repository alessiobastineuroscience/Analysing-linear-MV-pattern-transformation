function [Ttilde,optlambda,gof]=ridgeregmethod(x,y,lambdas); 
% Computation of the estimated linear transformation T (i.e. such that y=Tx) 
% via a cross-validated version of the ridge regression method.
% 
% INPUT
% x:         patterns in the ROIX
% y:         patterns in the ROIY
% lambdas:   set of possible regulariation parameter
% OUTPUT
% Ttilde:    estimated transformation
% optlambda: optimal regularization parameter
% gof:       goodness-of-fit
% Alessio Basti 20/02/2019 (Basti et al. 2019)

k=1;
for i=lambdas                    
   H=x'*pinv(x*x'+i*eye(size(x,1)))*x;                                      
   for j=1:size(x,2)
       A(j,j)=1/(1-H(j,j));
   end            
   CrossValid(k)=(norm(A*((eye(size(x,2))-H)*y'),'fro'))^2;                    
   k=k+1;
end
[B C]=min(CrossValid);
optlambda=lambdas(C);
gof=100*(1-B/(size(x,2)*size(y,1)));
Ttilde=y*x'*pinv(x*x'+optlambda*eye(size(x,1)));


return