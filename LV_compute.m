function [composition] = LV_compute(r_mean, S, D, A, T, step)


N_mean=1/(S*A*sqrt(pi/2)); 
sigma3=N_mean/12;
NN=zeros(T,S); 
r=ones(S,1);

% Interaction matrix
%   AA=exprnd(A,S,S);
    AA=rand(S,S)*A*2;
    AA=AA-diag(diag(AA))+eye(S); 

 % considering carrying capacity   
%  K=normrnd(1,sqrt(1/12),[S,1]);
%  for i=1:S
%      for j=1:S
%      AA(i,j)=AA(i,j)*K(j)/K(i);
%      end
%  end

 % Initial composition
  N0=abs(normrnd(N_mean,sigma3,1,S)); 
  NN(1,:)=N0;

 % Migration from source
   Migra=ones(1,S)*D;
 
% compute the dynamics in each patch
for i=2:T
       
       for j=1:S
           k1=r(j)*NN(i-1,j)*(1-AA(j,:)*(NN(i-1,:)'))*step;
           k2=r(j)*(NN(i-1,j)+k1/2)*(1-AA(j,:)*(NN(i-1,:)'))*step;
           k3=r(j)*(NN(i-1,j)+k2/2)*(1-AA(j,:)*(NN(i-1,:)'))*step;
           k4=r(j)*(NN(i-1,j)+k3)*(1-AA(j,:)*(NN(i-1,:)'))*step;
         NN(i,j)=NN(i-1,j)+(1/6)*(k1+2*k2+2*k3+k4);
      
       end

         NN(i,:)=NN(i,:)+Migra;
         
end

composition=NN;

end