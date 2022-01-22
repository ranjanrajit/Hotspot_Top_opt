function [senstivity_field] = senstivity_analysis_function(constant_term,T_fd,p,slab,use_xPhys,Temp_field_sn,Slab_thickness,no_of_slabs,con_penal,flux_penal,conductivity,Q_not,q_min,iKt,jKt,KE_thermal,L_one_slab)
%   This function takes adjoint information and finds senstivity
%% Gathering size data from K matrix and finding flux nodes
nelx=size(use_xPhys,2);
nely=size(use_xPhys,1);
local_node_nrs=(1+nelx)*(1+nely);
%% Adjoint equation intializations
lambda = zeros((nely+1)*(nelx+1),1);
%% Defining conductivity G matrix
K_min=conductivity*1e-4;
sKt = reshape(KE_thermal(:)*(K_min+use_xPhys(:)'.^con_penal*(conductivity-K_min)),16*nelx*nely,1);
G = sparse(iKt,jKt,sKt); G = (G+G')/2;
%% Defining Adjoint Load                                                                                                                                                                                                     
nodes_per_slab=(Slab_thickness+1)*(nelx+1);
Big_L=spalloc(no_of_slabs*(nelx+1),nodes_per_slab,(nelx+1));
Big_L(((slab-1)*(nelx+1))+1:(slab*(nelx+1)),:)=L_one_slab;
local_load=(T_fd.^(p-1))'*Big_L;
adj_load=(constant_term)*local_load'; % adjoint load
%% SOLVING
alldofs     = [1:(nely+1)*(nelx+1)];
fixeddofs=[1:nely+1:((nely+1)*nelx)+1];
freedofs    = setdiff(alldofs,fixeddofs);
lambda(freedofs,:) = G(freedofs,freedofs) \ adj_load(freedofs,:);   % finding lambdas as per Eq. (20) in journal paper
%% Define input heat flux and fixed bottom
% Flux scaled to conductivity
flux_ele=q_min+(Q_not-q_min)*(use_xPhys(nely,:).^flux_penal);
flux_ele=[0 flux_ele 0];
temp=movmean(flux_ele,2);
flux_node=temp(2:end);
q_in=[(nely+1):nely+1:((nely+1)*(nelx+1))];
Q(q_in)=flux_node;
%% Senstivity calculation as per equaion 19 in journal paper
temp_mat=zeros(nely+1,nelx+1);
temp_mat(nely+1,:)=((Q_not-q_min)*flux_penal)/2;
dQ=reshape(temp_mat,local_node_nrs,1);
const_mat=(conductivity-K_min)*con_penal*KE_thermal;
N_c=(Q_not*Slab_thickness)/conductivity;
T_s=Temp_field_sn*N_c;
for ely = 1:nely
    for elx = 1:nelx
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx   +ely;
        edof2 = [n1; n2; n2+1; n1+1];
        term1=(use_xPhys(ely,elx).^(con_penal-1))*const_mat*T_s(edof2);
        term2=((use_xPhys(ely,elx).^(flux_penal-1))*dQ(edof2))-term1;
        senstivity_field(ely,elx)=(lambda(edof2)')*term2;
    end
end
end

