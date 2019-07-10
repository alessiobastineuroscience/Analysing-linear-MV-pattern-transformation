function [rdd,density]=sparsityfeatures(Ttilde,x,y);     
% Computation of the density curve and of its rate of decay
% INPUT
% Ttilde:  estimated transformation
% x:       patterns in the ROIX
% y:       patterns in the ROIY
% OUTPUT
% rdd:     rate of decay of an exponential function fitted to the density curve
% density: density curve
% Alessio Basti 20/02/2019 (Basti et al. 2019)

Treshaped=reshape(Ttilde,1,prod(size(Ttilde)));             
for threshold=0:99
   Tthresholded=zeros(1,length(Treshaped));
   vec=find(abs(Treshaped)/(max(max(abs(Treshaped))))>(threshold/100));
   Tthresholded(vec)=Treshaped(vec);
   density(threshold+1)=length(vec)/length(Treshaped);
   Tthresholded=reshape(Tthresholded,size(Ttilde));
   var_noise(threshold+1)=norm(y-Tthresholded*x,'fro')/norm(y,'fro');    
end
func=fit((0:99)',100*density','exp1');
rdd=func.b;

return