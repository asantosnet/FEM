function [ test1xlast,test1ylast,test1x,test1y,test1xelastic,test1yelastic ] = HP_Modified_Newton(  )


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
delta = 0.006;
poisson = 0.3;
sigyield = 220e6;
t = 0.5;

nsteps = 300; % number of load steps

G = E/(2*(1+poisson));
K = E/(3*(1-2*poisson));

% Here we plot the mesh

[Edof,Ex,Ey,B1,~,B3,~,~,~,~,P4]=ringmesh(Ri,Ro,elemtype,hr,hphi);

% number of elements

nelement = size (Ex,1);

% number of degrees of freedom

ndofs = max(max(Edof));

AoF = zeros (ndofs,1); % Initial guess for a displacement

TOL = 10e-6; % Tolerance


% Dummy D

D = ones (4);

%  find the dofs for the constrained noeuds

cdofs = [B1(:,2);B3(:,2);P4(1,1)];

%  find the free dofs

fdofs = 1:ndofs;

fdofs(cdofs) = [];

% Identity matrix in voigt notation

Ivoigt = [1;
    1;
    1;
    0;
    0;
    0];


% 3x3 Identity matrix - voigt notation

Eyevoigtdev = [2/3 -1/3 -1/3 0 0 0;
    -1/3 2/3 -1/3 0 0 0;
    -1/3 -1/3 2/3 0 0 0;
    0 0 0 1/2 0 0;
    0 0 0 0 1/2 0;
    0 0 0 0 0 1/2;];

EyevoigtdevC = [2/3 -1/3 -1/3 0 0 0;
    -1/3 2/3 -1/3 0 0 0;
    -1/3 -1/3 2/3 0 0 0;
    0 0 0 2 0 0;
    0 0 0 0 2 0;
    0 0 0 0 0 2];



KS = zeros(ndofs,ndofs);
KTJ = zeros(ndofs,ndofs);

% The reaction forces

V = zeros(nsteps,1);

% TO save the values of g
% I won't... probably, make plus than 500 iterations per step...
gmatrix = zeros(500,nsteps);


AoF_1 = zeros (ndofs,1);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation - Modified Newton Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for step = 1:nsteps
    
    i = 1;
    
    g = ones (size(fdofs,1),1) ;
    
    % Input the boundary conditions
    
    AoF(B1(:,2),1) = -(delta/nsteps*step)/2;
    
    AoF(B3(:,2),1) = 0;
    
    AoF(P4(1,1),1) = 0;
    
    KTJ = zeros(ndofs,ndofs);
    
    % To obtain a better guess
    
    if step > 1
        
        AoF(fdofs) = 2*AoF(fdofs)  - AoF_1(fdofs);
        
    end
    
    AoF_1 = AoF;
    
    
    while norm(g) > TOL
        
        
        
        KS = zeros(ndofs,ndofs);
        
        EdF = extract(Edof,AoF);
        
        sigVM = zeros(nelement,1);
        
        for l = 1:nelement
            
            % the vector with all dofs associated with elements
            elemdofs = Edof(l,2:end);
            
            ep = [2,t];
            
            eq = [0;0];
            
            % retrieving es and et associated to each eple
            
            [~,et]=plants(Ex(l,:),Ey(l,:),ep,D,EdF(l,:));
            
            % Deviatoric strain - We use et as it has to be in voigt notation
            % Change et to six dim,
            
            etvoigt = [et(:,1);
                et(:,2);
                et(:,3);
                et(:,4);
                0;
                0];
            
            Etadeviatoric = Eyevoigtdev*etvoigt ;
            
            % the Von mises
            
            Etamises = sqrt(2/3)*sqrt(etvoigt'*Eyevoigtdev*etvoigt);
            
            % We will now define 1- G we will be using
            %                    2- Continuum tangent stifness Dt
            
            
            etacondition = sigyield/(3*G);
            
            if etacondition < Etamises
                Gstar =  sigyield/(3*Etamises);
                
                % The Dc
                
                Dt = 2*Gstar*Eyevoigtdev - Etadeviatoric*...
                    (4*sigyield/(9*(Etamises^3)))*Etadeviatoric'...
                    + K*(Ivoigt*Ivoigt');
            else
                Gstar = G;
                
                % The Dc
                
                Dt = 2*Gstar*Eyevoigtdev + K*(Ivoigt*Ivoigt');
            end
            
            % We will define here the Ds
            
            Ds = 2*Gstar*Eyevoigtdev + K*(Ivoigt*Ivoigt');
            
            % We will now calculate Ks
            
            [Ks,~]=plante(Ex(l,:),Ey(l,:),ep,Ds,eq);
            
            KS(elemdofs,elemdofs) = KS(elemdofs,elemdofs) + Ks;
            
            % We will now compute sigma
            
            sigma = Ds*etvoigt;
            
            sigVM(Edof(l,1)) = sqrt(3/2)*sqrt(sigma'*EyevoigtdevC*sigma);
            
            % we now compute the Jacobian
            
            switch i
                case 1
                    [KtJ,~]=plante(Ex(l,:),Ey(l,:),ep,Dt,eq);
                    KTJ(elemdofs,elemdofs) = KTJ(elemdofs,elemdofs) + KtJ;
                    % there is a initialization problem, the first iteration is
                    % not enough because of it. It convernges if we calculate KTJ up to
                    % the second iteration
                case 2
                    [KtJ,~]=plante(Ex(l,:),Ey(l,:),ep,Dt,eq);
                    KTJ(elemdofs,elemdofs) = KTJ(elemdofs,elemdofs) + KtJ;
                otherwise
            end
            
            
        end
        
        % solving the g(af) and we assume Ff to be zero
        
        Fint = KS*AoF;
        
        
        g = Fint(fdofs);
        gmatrix(i,step) = norm(g);
        
        %In order to update the Jacobian for the BFGS
        
        
        if norm(g) < TOL
            
            Q = KS*AoF;
            
            V(step) = -sum(Q(B1(:,2)));
            
            % if we want to know where we stopped to add values to gmatrix
            %( since the converging value migh be so small that we can
            % count on properly differentiating it from 0)
            % by using the sum we can be sure to know that there is no
            % conflict and that it will be the maximum value, i.e. we won't
            % by an accident choose a value that is
            % already there
            
            gmatrix (i+1,step) = sum(gmatrix(1:i,step));
            break
            
        end
        
        AoF(fdofs) = AoF(fdofs) - KTJ(fdofs,fdofs)\g;
        
        norm(g)
        
        i = i+ 1
        
        
    end
    
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure (8)


fill (Ex',Ey',sigVM');
shading flat
colorbar
colormap(jet);

figure(9)

plot(linspace((delta/(2*nsteps)),delta/2,nsteps),V(:,1),'b+');


xlabel('delta/2');
ylabel('Reaction Forces');

figure (10)

% we define the x and y axis, y being g(k+1) and x being g(k)


% We find the one who has the biggest amount of points

posMAXtest1x = 0;
test1x = 0;
stepmax = 0;
for gamma = 1:nsteps
    test = gmatrix (:,gamma);
    posMAXx =find(test == max(test));
    
    if posMAXx > posMAXtest1x
        posMAXtest1x = posMAXx;
        test1x = test;
        stepmax = gamma;
    end
    
end



% we use the previouus vector to define the vector on the y axis
test1y = test1x(2:(posMAXtest1x-1));

% we remove the unwanted values from test1x
test1x((posMAXtest1x-1):end) = [];

% we will plot for the last step
test1xelastic = gmatrix (:,1);

% we find the number of iterations that were used
posMAXtest1xelastic =find(test1xelastic == max(test1xelastic));

% we use the previouus vector to define the vector on the y axis
test1yelastic = test1xelastic(2:(posMAXtest1xelastic-1));

% we remove the unwanted values from test1x
test1xelastic((posMAXtest1xelastic-1):end) = [];

% we will plot for the last step
test1xlast = gmatrix (:,nsteps);

% we find the number of iterations that were used
posMAXtest1xlast =find(test1xlast == max(test1xlast));

% we use the previouus vector to define the vector on the y axis
test1ylast = test1xlast(2:(posMAXtest1xlast-1));

% we remove the unwanted values from test1x
test1xlast((posMAXtest1xlast-1):end) = [];

% Setting the axis as equals - gca is the current axis handle



loglog(test1x,test1y,'-g');

xLimits = [TOL*10e-3 1];
yLimits = [TOL*10e-3 1];
axes_equal_loglog (gca,xLimits,yLimits)

xlabel('|g(k)|');
ylabel('|g(k+1)|');

hold on

% Setting the axis as equals - gca is the current axis handle

loglog(test1xlast,test1ylast,'-s');
loglog(test1xelastic,test1yelastic,'-r');

grid on
legend(strcat('For step = ',num2str(stepmax)),strcat(' For last step = ',...
    num2str(nsteps)),'Elastic domain');

hold off




end

