function [T_hat] = hotspot_analysis_function(use_xPhys,con_penal,flux_penal,conductivity,Q_not,q_min,iKt,jKt,KE_thermal)
%   This function takes slab density as input
%   and performs FEA assuming a penalized conductivity heat flux from top 
%   and bottom nodes as sink. This fuction return the normalized temperature 
%   field as output using thr normalization constamt as N_c=(q_0*s)/(k_0)
%% Gathering size data from density field
nelx=size(use_xPhys,2);
nely=size(use_xPhys,1);
K_min=conductivity*1e-4;
%% Preparing the conductivity matrix (G)
Q = sparse((nely+1)*(nelx+1),1); 
T = zeros((nely+1)*(nelx+1),1);
sKt = reshape(KE_thermal(:)*(K_min+use_xPhys(:)'.^con_penal*(conductivity-K_min)),16*nelx*nely,1);
G = sparse(iKt,jKt,sKt); G = (G+G')/2;
%% Define input heat flux and fixed bottom
% Flux scaled to conductivity
flux_ele=q_min+((Q_not-q_min)*(use_xPhys(nely,:).^flux_penal));
flux_ele=[0 flux_ele 0];
temp=movmean(flux_ele,2);
flux_node=temp(2:end);
q_in=(nely+1):nely+1:((nely+1)*(nelx+1));
Q(q_in)=flux_node;
% fixed bottom
fixeddofs=[1:nely+1:((nely+1)*nelx)+1];
alldofs     = [1:(nely+1)*(nelx+1)];
freedofs    = setdiff(alldofs,fixeddofs);
%% SOLVING
T(freedofs,:) = G(freedofs,freedofs) \ Q(freedofs,:);   
T(fixeddofs,:)= 0;
N_c=(Q_not*nely)/conductivity;
T_hat=T./N_c;    % Normalized temperature array
end

