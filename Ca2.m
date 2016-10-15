close all
clear all
clc


% Defining constants
% Ri    - Inner radius of the ring
% Ro    - Outer radius of the ring
% elemtype - Element type    "tria3"  - 3-node triangular element
%                            "quad4" - 4-node quadrilateral element
% hr    - Approximate element size in radial direction
% hphi  - Approximate element size in circumferencial direction

Ri  = 0.15;
Ro = 0.25 ;
elemtype = 'tria3';
hr = 0.1;
hphi = 0.1;
E = 210e9 ;
ndelta = 100;
maxdelta = 0.008;
mindelta = 0.001;
poisson = 0.3;
sigyield = 220e6;
t = 0.5;


% Here we plot the mesh

[Edof,Ex,Ey,B1,B2,B3,B4,P1,P2,P3,P4]=ringmesh(Ri,Ro,elemtype,hr,hphi);

% number of elements

nelement = size (Ex,1);

% type 1   - axysimmetric

D = hooke(2,E,poisson);


% Calculating the boundary Conditions
% The boundary condiions for B1 - moves - we consider here y not moving,
% symmetry, thats why we take only the values for y degrees of freedom
% now we do the same thing for the other side
% The boundary condiions for B3 - doesn't move - y doesnt move
% we gotta also say that the x doesn't moves, meaning
% calculate K for all elements


deltavect = linspace(mindelta,maxdelta,ndelta); % Displacement vector

% We do all that until we have max stress above yield stress

gamma = 1;

sigmax = 0;

% Vertical reaction

V = zeros (ndelta,1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while (sigyield + 100e6) > sigmax
    
    
    delta = deltavect (gamma);
    
    % find the dofs for the constrained noeuds
    
    cdofs = [B1(:,2);B3(:,2);P4(1,1)];
    
    % set the bc's
    
    bc = zeros(length(cdofs),2);
    bc(:,1) = cdofs;
    bc(1:size(B1,1),2) = -delta/2 ;
    
    % find the free dofs
    
    ndofs = max(max(Edof));
    
    fdofs = 1:ndofs;
    
    fdofs(cdofs) = [];
    
    % We remove values existing in dofs that match the ones existing in cdofs
    
    K = zeros (ndofs,ndofs);
    F = zeros (ndofs,1);
    
    % we will here do the partitioning of K
    
    for n = 1:nelement
        
        % the vector with all dofs associated with elements
        
        elemdofs = Edof(n,2:end);
        
        ep = [1,t];
        eq = [0;0];
        
        % retrieving K and Fe associated with each element
        
        [ke,fe]=plante(Ex(n,:),Ey(n,:),ep,D,eq);
        
        % Assembling everything
        
        K(elemdofs,elemdofs) = K(elemdofs,elemdofs) + ke;
        F(elemdofs) = F(elemdofs) + fe;
        
        
    end
    
    
    % Here we recover Kff
    
    Kff = K(fdofs,fdofs);
    
    % Here we recover Kfc
    
    Kfc = K(fdofs,cdofs);
    
    % Here we recover kcf
    
    Kcf = K(cdofs,fdofs);
    
    % here we define Kcc
    
    Kcc = K(cdofs,cdofs);
    
    A = zeros(ndofs,1);
    
    A(B1(:,2)) = -delta/2;
    
    % Here we calculate the displacement for the free node
    
    A(fdofs) = (Kff)\(F(fdofs) - Kfc*bc(:,2));
    Q=K*A;
    Ed = extract (Edof,A);
    
    
    % here we calculate the element stress  matrix and element strain matrix
    
    sigVM = zeros(nelement,1);
    sigmaxvect = zeros(nelement,1);
    
    for l = 1:nelement
        
        % the vector with all dofs associated with elements
        
        elemdofs = Edof(l,2:end);
        ep = [1,t];
        eq = [0;0];
        
        % retrieving es and et associated to each eple
        
        [es,et]=plants(Ex(l,:),Ey(l,:),ep,D,Ed(l,:));
        
        % our stress matrix will be
        
        S = [es(:,1) es(:,4) 0;
            es(:,4) es(:,2) 0;
            0 0 es(:,3)];
        
        sigmaxvect(l) = max(max(S));
        
        % here we define the mean stress
        meanstress = (1/3)*trace(S);
        
        % we calculate the Deviatoric stress
        
        sigdev = S - meanstress*eye(3,3);
        
        % we define the effective stress - Von Mises
        
        sigVM(Edof(l,1)) = sqrt(3/2)*norm(sigdev);
        
        ny = [0;1;0];
        
        
        
        
        
        
    end
    
    % we calculate here the vertical reaction
    
    ty = Q(B1(:,2));
    % to get vertical reaction force
    V(gamma) = sum(ty);
    sigmax = max(sigmaxvect);
    
    
    gamma = gamma + 1;
    
    
end

% Removing if zero

V(V==0)=[];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)


plotpar = [1 4 1];
eldisp2(Ex,Ey,Ed,plotpar,100)



figure (2)


fill (Ex',Ey',sigVM');
shading flat
colorbar
colormap(jet);

figure (3)

plot(deltavect(1:size(V,1))/2,V(:,1),'-b');


xlabel('delta/2');
ylabel('Reaction Forces');

