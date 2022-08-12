clear;clc;close all;
r_mean=1.0; % mean growth rate
S=1:1:100; % number of species 
DD=1e-6; % migration rate
A_mean_range=[0.01 1.5]; % mean inhibition alpha
ystep=100; %for A_mean
T_real=1000; %total evolution time
step=0.05; % time interval 
T=ceil(T_real/step); %number of circulation
num_sim=100; % number of simulation
period=200; % duration for calculating richness
% composition_full=zeros(length(S)*ystep*num_sim,length(S)); % record the final compositions of all simulations
abundance_threshold=1e-3; %threshold for survival species

fluc_record=zeros(length(S)*ystep,num_sim);
diversity_record=zeros(length(S)*ystep,num_sim); 
A_mean=linspace(A_mean_range(1),A_mean_range(2),ystep); 

for mm=1:length(S)
    for pp=1:ystep
      for hh=1:num_sim
         composition= LV_compute(r_mean, S(mm), DD, A_mean(pp),T,step);  
         diversity_index=max(composition(T-period+1:T,:)); richness=find(diversity_index>abundance_threshold);
         diversity_record((mm-1)*ystep+pp,hh)=length(richness)/S(mm);
         fluc_CV=std(composition(T-500:T,:))./mean(composition(T-500:T,:)); fluc_CV(isnan(fluc_CV))=0;
         fluc_record((mm-1)*ystep+pp,hh)=mean(fluc_CV(fluc_CV>0));
      end
    end
    
end


