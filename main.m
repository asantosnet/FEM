close all
clear all


% Initialization of constants

rhowater = 1000;  %%
g = 9.81;  %% m/s^2
rhoconcrete = 2300;  %%
poisson = 0.18;
E = 22e9;  %%
nheight = 20  ;  %% m
nthick = 10;   %% m
waterheight = 200;  %% m
Patm = 1e5 ; %%


% type  3 - axysimmetric

D = hooke(3,E,poisson);

% recovering positions and normal values

[Ex,Ey,Edof,Boundarygeom,Boundarydof]=hoover2d(nthick,nheight);

% number of boundary segments

nbound = size (Boundarydof);

% number of degrees of freedom

ndofs = max(max(Edof(:, 2:end)));


f_b = zeros(ndofs,1);

% matrix N that will be used later to calculate f_b for the boundary where
% we have water preassure

N = [ 1,0;
    0,1;
    1,0;
    0,1;
    ];

% calculating the stresses/ boundary conditions for water side. As for the
% air side and goround we consider the preassure to be vector T to be zero
% we don't need to calculate them, as it was already done when we
% initialized f_b

for n = 1:nbound(1)
    
    x1 = Boundarygeom(n,1);
    x2 = Boundarygeom(n,3);
    y1 = Boundarygeom(n,2);
    y2 = Boundarygeom(n,4);
    nx = Boundarygeom(n,5);
    ny = Boundarygeom(n,6);
    
    lf =sqrt((x1-x2)^2+(y1-y2)^2);
    
    Normal = [nx;ny];
    
    dofs = Boundarydof(n,:);
    fbf = zeros(4,1);
    
    
    % for water side
    
    if (ny ~= 1) &&(ny ~= -1) && (y1 || y2 < waterheight) && (sign(nx)==1)
        
        % 1 by 2 matrix with the traction on the water side, it will be
        % zero otherwise
        
        T  = -((rhowater*g*(waterheight-(y1+y2)/2)))*(Normal);
        fbf = (lf*(N)*(x1+x2)*T)/4;
        
    end
    
    f_b(dofs) = f_b(dofs) + fbf(:);
    
end

% boundary conditions for the botton

bc = zeros (nthick*2+2,2);
bc (:,1) = 1:(nthick*2 +2);

% number of elements

nelement = size (Ex);

K = zeros (ndofs,ndofs);
F = zeros (ndofs,1);

% Calculating the K and F matrix for each element

for k = 1: nelement(1)
    
    
    ex = Ex(k,:);
    ey = Ey(k,:);
    eq = [0;0];
    
    % the bodyforce
    
    eq(2) = -rhoconcrete*g;
    
    % n x n integration points where n the number of gauss poits
    % xa and ya = 0 ,nxa = 0 and nya = 1
    
    ep = [0;0;0;1;2];
    
    [ke,fe]=axi4e(ex,ey,ep,D,eq);
    
    [K,F] = assem(Edof(k,:),K,ke,F,fe);
    
end

% Calculating the displacements for each element

F = F + f_b;

[Af] = solveq(K,F,bc);

% This will be usde to plot the unmodified dam
sizeF = size(F);
Fzero = zeros (sizeF(1),sizeF(2));

[Afzero] = solveq(K,Fzero,bc);

Edzero = extract (Edof,Afzero);

% Extract the displacement for each element

Ed = extract (Edof,Af);

Sprincipal = zeros(3,nelement(1));

% Calculating the principal stresses in each element

for l = 1: nelement(1)
    
    
    % xa and ya = 0
    
    ep = [0;0;0;1;2];
    
    
    % here we obtain a matrix with the stress associated with each point
    
    [es,~,~]=axi4s(Ex(l,:),Ey(l,:),ep,D,Ed(l,:));
    
    % here we calculate the general matrix by taking the mean of all gauss
    %
    sprincipal = [mean(es(:,1)),0,mean(es(:,4));
        0,mean(es(:,3)),0;
        mean(es(:,4)),0,mean(es(:,2));
        ];
    
    
    Sprincipal(:,l) = eig(sprincipal);
    
end

% maximum stress    ;

max_Sprincipal  = zeros (nelement(1),4);

max_Sprincipal(:,1) = (Sprincipal(3,:));

% minimum stress

min_Sprincipal  = zeros (nelement(1),4);

min_Sprincipal(:,1) = (Sprincipal(1,:));


figure (1)


fill (Ex',Ey',min_Sprincipal');
shading flat
colorbar
colormap(jet);
figure (2)

fill (Ex',Ey',max_Sprincipal');
shading flat
colorbar
colormap (jet)

figure(3)

hold on

plotpar = [1 4 1];
eldisp2(Ex,Ey,Ed,plotpar,1000)

plotpar = [1 2 1];
eldisp2(Ex,Ey,Edzero,plotpar,1000)

hold off
