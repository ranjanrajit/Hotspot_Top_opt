function [ Temp_field_domain,sens_domain,norm_p] = hotspot_analysis(xPhys,T_critical,p,no_of_nodes_overlap,Slab_thickness,con_penal,flux_penal,conductivity,Q_not,q_min,factor,iKt,jKt,KE_thermal,L_one_slab)
%% Adding the substrate plate at bottom
nelx=size(xPhys,2);
nely=size(xPhys,1);
substrate_plate_down=ones(Slab_thickness,nelx);
xPhys=[substrate_plate_down;xPhys];
no_of_slabs=nely+1;
%% Doing the hotspot analysis for every slab in parallel
parfor slab=1:1:no_of_slabs
    use_xPhys=xPhys(slab:slab+Slab_thickness-1,:);  %this one is elemental
    [Temp_field_sn(:,slab)]=hotspot_analysis_function(use_xPhys,con_penal,flux_penal,conductivity,Q_not,q_min,iKt,jKt,KE_thermal);
end
%% Doing the senstivity analysis for each slab 
top_nodes = (Slab_thickness+1):Slab_thickness+1:((Slab_thickness+1)*(nelx+1));
all_top_temp = Temp_field_sn(top_nodes,:);
T_fd=all_top_temp(:);
g=T_fd./T_critical;
g_p=g.^p;
sig=sum(g_p);
norm_p=(sig/no_of_nodes_overlap)^(1/p);
const_term_num=factor*(sig^((1/p)-1));
denom1=no_of_nodes_overlap^(1/p);
denom2=T_critical^p;
N_c=(Q_not*Slab_thickness)/conductivity;
const_term_denom=denom1*denom2*N_c;
constant_term=const_term_num/const_term_denom;
parfor slab=1:1:no_of_slabs
    use_xPhys=xPhys(slab:slab+Slab_thickness-1,:);  %this one is elemental
    senstivity_field{slab}=senstivity_analysis_function(constant_term,T_fd,p,slab,use_xPhys,Temp_field_sn(:,slab),Slab_thickness,no_of_slabs,con_penal,flux_penal,conductivity,Q_not,q_min,iKt,jKt,KE_thermal,L_one_slab);
end
%% Assembling all the slabs and senstivities
for slab=1:1:no_of_slabs
    global_temp_field(slab:slab+Slab_thickness,:,slab)=reshape(Temp_field_sn(:,slab),Slab_thickness+1,nelx+1); %nodal
    global_senstivity_field(slab:slab+Slab_thickness-1,:,slab)=senstivity_field{slab}; %elemental
end
df=sum(global_senstivity_field,3);
%% Getting temp field for visualization purpose
Temp_field_max=max(global_temp_field,[],3);
%% Passing back only the domain data (i.e. excluding substrate)
Temp_field_domain=Temp_field_max(Slab_thickness+1:nely+Slab_thickness+1,:);
sens_domain=df(Slab_thickness+1:nely+Slab_thickness,:);
end
