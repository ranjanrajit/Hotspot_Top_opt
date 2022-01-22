%% Code Description
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves the standard compliance minimization problem           %
% for a cantilever loading case with the hotspot constraint to avoid local%
% overheating in AM process. The code is written by Rajit Ranjan, PhD     %
% candidate in The Department of Precision and Microsystems Engineering   %
% (PME) at TU Delft, Netherlands. Please send your comments to the author %
% at: r.ranjan@tudelft.nl                                                 %
%                                                                         %
% A novel constraint with hotspot avoidance is added. The code is based   %
% on the 88 line top opt code provided by Andreassen et al. (2011).       % 
% In addition to this, the robust formulation by Wang et al. (2011) is    %
% used.                                                                   %
%                                                                         %
% Three new functions are added which perform the hotspot analysis        %
% and calculate sensitivities. The names of these functions are           %
% as follows:                                                             %
% 1. hotspot_analysis.m                                                   %
% 2. hotspot_analysis_function.m                                          %
% 3. senstivity_analysis_function.m                                       %
%                                                                         %
% The section of code which is newly                                      % 
% added for hotspot TO is labeled as 'NEW FOR HOTSPOT TO'                 %
%                                                                         %
% The novel part of the code is commented for users                       %
% to follow and variable names are mostly consistent with the journal     %
% paper. The code can work with method of moving asymptotes (MMA) for     %
% optimization and users are encourged to plugin their own MMA code.      %
% Comments are provided where MMA code needs to be added. Due to this,    %
% the code cannot be directly used.                                       %
%                                                                         %
% Disclaimer:                                                             %
% The author reserves all rights but does not guaranty that the code is   %
% free from errors. Furthermore, he shall not be liable in any event      %
% caused by the use of the program.                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CLEAR THE WORKSPACE
clc
clearvars
%% INITIALIZE
nelx=100;               % no of element in the X axis
nely=60;                % no of element in the Y axis
volfrac=0.5;
penal=3;
rmin=6;
ft=2;
p=15;                   % P value for defining P-norm
%% MATERIAL PROPERTIES (STRUCTURAL ANALYSIS)
E0 = 1;
Emin = 1e-9;
nu = 0.3;
%% MATERIAL PROPERTIES (THERMAL ANALYSIS) (NEW FOR HOTSPOT TO)
conductivity=1;         % represented as k_0 in journal paper
%% SLAB ANALYSIS PARAMETERS (NEW FOR HOTSPOT TO)
Slab_thickness=12;      % provide slab thickness s here in terms of number of elements
T_critical=2.1;         % provide T_critical value here based on calibration, 2.1 for 45 deg
Q_not=1;                % input flux
K_min=1e-4;             % mimimum value of conductivity
q_min=0;                % minimum value of heat flux
con_penal=3;            % thermal penalization exponenet, represented as r in the journal paper
flux_penal=3;             % thermal penalization exponenet, represented as r in the journal paper
factor=1;               % for p_norm scaling based on the continuation scheme
%% PREPARE FINITE ELEMENT ANALYSIS (STRUCTUTAL)
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE_structural = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%% PREPARE FINITE ELEMENT ANALYSIS (THERMAL) (NEW FOR HOTSPOT TO)
KE_thermal = [ 2/3 -1/6 -1/3 -1/6
      -1/6  2/3 -1/6 -1/3
      -1/3 -1/6  2/3 -1/6
      -1/6 -1/3 -1/6 2/3];
row=1;
for elx = 1:nelx
    for ely = 1:Slab_thickness
        n1 = (Slab_thickness+1)*(elx-1)+ely;
        n2 = (Slab_thickness+1)* elx   +ely;
        edofMat_th_slab(row,:) = [n1+1; n2+1; n2; n1];
        row=row+1;
    end
end
iKt = reshape(kron(edofMat_th_slab,ones(4,1))',16*nelx*Slab_thickness,1);
jKt = reshape(kron(edofMat_th_slab,ones(1,4))',16*nelx*Slab_thickness,1);
%% DEFINE LOADS AND SUPPORTS (CANTILEVER)
F = sparse(2*(nely+1)*(nelx+1),1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
xmin = 1e-6;
n = nelx*nely; m = 2;
x = repmat(volfrac,nely,nelx);

xTilde = x;
% ROBUST PARAMETERS
beta = 2;
etad = 0.25;
etai = 0.5;
etae = 1-etad;
xPhysd = (tanh(beta*etad) + tanh(beta*(xTilde - etad)))./(tanh(beta*etad) + tanh(beta*(1 - etad)));
xPhysi = (tanh(beta*etai) + tanh(beta*(xTilde - etai)))./(tanh(beta*etai) + tanh(beta*(1 - etai)));
xPhyse = (tanh(beta*etae) + tanh(beta*(xTilde - etae)))./(tanh(beta*etae) + tanh(beta*(1 - etae)));

vd = mean(xPhysd(:));
vi = mean(xPhysi(:));
ve  =  mean(xPhyse(:));
vds = vd;

xminvec  = xmin*ones(n,1);
xmaxvec  = ones(n,1);
low   = xminvec;
upp   = xmaxvec;
cmma = 10000*ones(m,1); % Belang van constraint tov objective often 1000, 10000
d = zeros(m,1);
a0 = 1;
a = zeros(m,1);
xold1 = x;
xold2 = x;
%% TOP NODES TEMPERATURE EXTRACTION MATRIX (NEW FOR HOTSPOT TO)
% only topmost nodes are taken for finding the maximum
% The matrix called 'L_one_slab' extracts top node temperatures
no_of_nodes_overlap=(Slab_thickness)*(nely+1.5-(Slab_thickness/2))*(nelx+1);
nodes_per_slab=(Slab_thickness+1)*(nelx+1);
top_nodes = (Slab_thickness+1):Slab_thickness+1:((Slab_thickness+1)*(nelx+1));
L_one_slab = spalloc((nelx+1),nodes_per_slab,nelx+1);
for i=1:1:nelx+1
    L_one_slab(i,top_nodes(i))=1;
end
%% LOOP STARTS
loop = 0;
change = 1;
beta_max = 128;
M_nd = 100;
counter=1;
%% START ITERATION
while loop<401
  loop = loop + 1;
  %% FE-ANALYSIS
  % only eroded design is considered for compliance evaluation
  sK = reshape(KE_structural(:)*(Emin+xPhyse(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
  K = sparse(iK,jK,sK); K = (K+K')/2;
  U(freedofs) = K(freedofs,freedofs)\F(freedofs);
  %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
  ce = reshape(sum((U(edofMat)*KE_structural).*U(edofMat),2),nely,nelx);
  c = sum(sum((Emin+xPhyse.^penal*(E0-Emin)).*ce));
  dce = -penal*(E0-Emin)*xPhyse.^(penal-1).*ce;
  dv = ones(nely,nelx);
  %% APPLYING P-NORM SCALING FACTOR (NEW FOR HOTSPOT TO)
  if rem(loop,25)==0
     max_g=max(Temp_field(:))/T_critical;
     factor=max_g/norm_p;
  end
  %% SLAB AND SENSTIVITY ANALYSIS (NEW FOR HOTSPOT TO)
  % only intermediate design is considered for hotspot analysis
  % The function 'slab_analysis' takes current intermediate design as input
  % and returns hotspot field (variable name: Temp_field), senstivity 
  % field (variable name: senstivity) and P-norm value (variable name: 
  % norm_p) as output
  [Temp_field,senstivity,norm_p]=hotspot_analysis(flipud(xPhysi),T_critical,p,no_of_nodes_overlap,Slab_thickness,con_penal,flux_penal,conductivity,Q_not,q_min,factor,iKt,jKt,KE_thermal,L_one_slab);
  senstivity=flipud(senstivity);
  %% THERMAL CONSTRAINT (NEW FOR HOTSPOT TO)
  tru_max=factor*norm_p;
  %% PROJECTION FILTER SENSTIVITY
   dxxd = (beta*(sech(beta*(xTilde-etad))).^2)./(tanh(beta*etad) + tanh(beta*(1-etad)));
   dxxi = (beta*(sech(beta*(xTilde-etai)).^2))./(tanh(beta*etai) + tanh(beta*(1-etai)));
   dxxe = (beta*((sech(beta*(xTilde-etae))).^2))./(tanh(beta*etae) + tanh(beta*(1-etae)));
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft == 1
    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
  elseif ft == 2
    dce(:) = H*((dce(:).*dxxe(:))./Hs);
    dv(:) = H*((dv(:).*dxxd(:))./Hs);
    senstivity(:) = H*((senstivity(:).*dxxi(:))./Hs);
  end
   

  f0val = c;
  df0dx = dce(:);
  df0dx2 = 0*df0dx;
  fval(1,:) = sum(xPhysd(:))/(n*vds) -1;
  fval(2,:) = tru_max-1;
  dfdx(1,:) = dv(:)'/(n*vds);
  dfdx(2,:) = senstivity(:)';                                                                                                                                   
  dfdx2 = 0*dfdx;
  xval = x(:);

  %% MMA DESIGN UPDATE:                    
  % ADD MMA (or any other optimizer) to get the updated design variables:
  % xnew

  xTilde(:) = (H*xnew(:))./Hs;
  
  xPhysd = (tanh(beta*etad) + tanh(beta*(xTilde - etad)))./(tanh(beta*etad) + tanh(beta*(1 - etad)));
  xPhysi = (tanh(beta*etai) + tanh(beta*(xTilde - etai)))./(tanh(beta*etai) + tanh(beta*(1 - etai)));
  xPhyse = (tanh(beta*etae) + tanh(beta*(xTilde - etae)))./(tanh(beta*etae) + tanh(beta*(1 - etae)));
  
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  vd = mean(xPhysd(:));
  if(mod(loop,20)==1)
      vds = volfrac/mean(xPhysi(:))*vd;
  end
  %% updating beta....
   if (mod(loop,50)==0 && beta<=beta_max)
       beta = 2*beta;
   end
   if(beta>beta_max)
       beta = beta_max;
   end
   %% SAVING DATA FROM EACH ITERATION (NEW FOR HOTSPOT TO)
    M_nd1=4.*xPhysi.*(1-xPhysi);
    M_nd=(sum(sum(sum(M_nd1)))/n)*100;
    cons_vec(loop) = fval(1,:);
    c_vec(loop) = c;
    M_nd_vec(loop) = M_nd;
    %% PLOTTING
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f %f\n',loop,c, ...
        mean(xPhysi(:)),change,tru_max);
    figure(1)
    colormap(gray); imagesc(1-xPhysi); caxis([0 1]); axis equal; axis off; drawnow;
    figure(2)
    colormap(jet); imagesc(flipud(Temp_field));colorbar; axis equal; axis off; drawnow; 
end
%% CONVERGANCE PLOT
vector = 1:1:loop-1;
plot(vector,c_vec,'Linewidth',2)
xlabel('Iterations')
ylabel('compliance')