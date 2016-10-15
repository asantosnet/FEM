%
%
%
%
%
%
close all
clear all

plate = 1; % Kirchoff

ptype = 1;
t = 0.008; % m
Lx = 3; % m
Ly = 1; % m
E = 210e9; % PA
v = 0.3;
qz =1000; % in Newton / m^ 2
tractiony = -5e6; % Newton / m^ 2
irb = 2;
irs = 1;
ep = [ptype t irb irs];
eq = [0 0 qz];

nx = 20;
lx = Lx/nx;
ny = 10;


D = hooke (ptype,E,v);

G = (E*eye(2,2))/(2*(1+v)); % 



[ Coord, Dof, Enode, Edof, Ex, Ey, DofRight, DofTop, ...
    DofLeft, DofBottom, cdof ] = quadmesh( Lx, Ly, nx, ny, 'no' );


nelement = size(Ex,1);
ndofs = max(max(Edof));

% find the dofs for the constrained noeuds
% For the Dofs along z where we have boundary conditions along Top and
% Bottom boundaries

DofBoundary = [DofBottom(:,1);
    DofTop(:,1)];


DofBoundaryZ = DofBoundary(3:3:end,1);

constraineddofs = [DofLeft(:,1);
    DofRight(:,1);
    DofBoundaryZ(:,1)];


% Find the free dofs

fdofs = 1:ndofs;
fdofs(constraineddofs) = [];

% set the bc's

bc = zeros(length(constraineddofs),2);
bc(:,1) = constraineddofs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate A and Ed -Deformation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[A,Ed,F]  = displacementcalc( ndofs,nelement,Edof, Ex, Ey,  DofTop, ...
    ep,eq,D,plate,G,bc,fdofs,constraineddofs,tractiony*ep(2),lx);
% PLotting the deformations
scx = 100;
scy = 200;
scz = 200;

xdofs = [1 6 11 16];
ydofs = [2 7 12 17];
zdofs = [3 8 13 18];

figure;
Dispviz3(Ex + Ed(:,xdofs)*scx,Ey + Ed(:,ydofs)*scy,Ed(:,zdofs)*scz,Edof,A,3,2);
colorbar
shading flat
colormap(jet);
title('3D plot of  the Deformation - Kirchoff');
xlabel( colorbar,'Deformation(m)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

es=zeros(nelement,3);
et=zeros(nelement,3);
ef=zeros(nelement,8);
ec=zeros(nelement,3);
s_min = zeros(nelement,1);
Mxx=zeros(nelement,1);
Myy=zeros(nelement,1);
Mxy=zeros(nelement,1);



for teta=1:nelement
    [es(teta,:),et(teta,:),ef(teta,:),ec(teta,:)]=shell2rs(Ex(teta,:)...
        ,Ey(teta,:),[ep(1) ep(2)],D,Ed(teta,:));
    
    
    S=[es(teta,1) es(teta,3) 0;    %Building of the stress matrix
        es(teta,3) es(teta,2) 0;
        0 0 0];
    
    eigen=eig(S);            %principal stresses
    
    s_min(teta)=min(eigen);
    Mxx(teta)=ef(teta,1);          %value of Mxx Myy Mxy
    Myy(teta)=ef(teta,2);
    Mxy(teta)=ef(teta,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second Part calculation using A and Ed / Buckling
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We need to use kirchoff
plate2nd = 1;

% Just so we recall what are we using
ep = [ptype t irb irs];
eq = [0 0 qz];


[ X,~,posvec,freedof,~] = bucklingcalc(Ex,Ey,ep,D,es,G,plate2nd...
                                       ,eq,bc,ndofs,Edof,nelement );
                                   
% Retrieving the data for each mode usign the posvec, where we will obtain 
% the position of each lambada in a crescenting order

u_mode = zeros(ndofs,6);
u_mode(freedof,:) = X(:,posvec(:,1));

Ed_mode1 = extract(Edof,u_mode(:,1));
Ed_mode2 = extract(Edof,u_mode(:,2));
Ed_mode3 = extract(Edof,u_mode(:,3));
Ed_mode4 = extract(Edof,u_mode(:,4));
Ed_mode5 = extract(Edof,u_mode(:,5));
Ed_mode6 = extract(Edof,u_mode(:,6));

% Plotting

scx = 50;
scy = 20;
scz = 20;

xdofs = [1 6 11 16];
ydofs = [2 7 12 17];
zdofs = [3 8 13 18];

figure;
Dispviz3(Ex + Ed_mode1(:,xdofs)*scx,Ey + Ed_mode1(:,ydofs)*scy,Ed_mode1(:,zdofs)*scz,Edof,u_mode(:,1),3,2);
colorbar
colormap(jet);
title('3D plot of  the Deformation - Buckling using Kirchoff mode 1');
xlabel( colorbar,'Deformation(m)');

figure;

Dispviz3(Ex + Ed_mode2(:,xdofs)*scx,Ey + Ed_mode2(:,ydofs)*scy,Ed_mode2(:,zdofs)*scz,Edof,u_mode(:,2),3,2);
colorbar
colormap(jet);
title('mode 2');
xlabel( colorbar,'Deformation(m)');

figure;

Dispviz3(Ex + Ed_mode3(:,xdofs)*scx,Ey + Ed_mode3(:,ydofs)*scy,Ed_mode3(:,zdofs)*scz,Edof,u_mode(:,3),3,2);
colorbar
colormap(jet);
title(' mode 3');
xlabel( colorbar,'Deformation(m)');

figure;

Dispviz3(Ex + Ed_mode4(:,xdofs)*scx,Ey + Ed_mode4(:,ydofs)*scy,Ed_mode4(:,zdofs)*scz,Edof,u_mode(:,4),3,2);
colorbar
colormap(jet);
title(' mode 4');
xlabel( colorbar,'Deformation(m)');

figure;

Dispviz3(Ex + Ed_mode5(:,xdofs)*scx,Ey + Ed_mode5(:,ydofs)*scy,Ed_mode5(:,zdofs)*scz,Edof,u_mode(:,5),3,2);
colorbar
colormap(jet);
title(' mode 5');
xlabel( colorbar,'Deformation(m)');

figure;

Dispviz3(Ex + Ed_mode6(:,xdofs)*scx,Ey + Ed_mode6(:,ydofs)*scy,Ed_mode6(:,zdofs)*scz,Edof,u_mode(:,6),3,2);
colorbar
colormap(jet);
title('6');
xlabel( colorbar,'Deformation(m)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Last part Kirchoff vs Mindlin plate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deflectionKtotal = zeros(50,1);
deflectionMtotal1 = zeros(50,1); % 2 and 1
deflectionMtotal2 = zeros(50,1); % 2 and 2
deflectionMtotal3 = zeros(50,1); % 2 and 3

counter = 1;
% In order to change the thickness, it will fo from 0.001 to 0.5001 in an
% 0.01 interval

for i= 0.001:0.01:0.5001

    ep = [ptype i 2 1];
    
    % calculate for Kirchoff

    [ AKK,EdKK ] = deflectioncalc( ndofs,nelement, Edof, Ex, Ey,  DofTop, ...
        ep,eq,D,1,G,bc,fdofs,constraineddofs,tractiony,cdof,lx );

    deflectionKtotal(counter,1) = (i^3)*sum(EdKK);

    % Calculate for Mindling and irb = 2 nad irs = 1

    [AMM1, EdMM1 ] = deflectioncalc( ndofs,nelement, Edof, Ex, Ey,  DofTop, ...
        ep,eq,D,2,G,bc,fdofs,constraineddofs,tractiony,cdof,lx);

     deflectionMtotal1(counter,1) = (i^3)*sum(EdMM1);

    % Calculate for Mindling and irb = 2 nad irs = 2
     
    ep = [ptype i 2 2];

    [ AMM2,EdMM2 ] = deflectioncalc( ndofs,nelement, Edof, Ex, Ey,  DofTop, ...
        ep,eq,D,2,G,bc,fdofs,constraineddofs,tractiony,cdof,lx);

    deflectionMtotal2(counter,1) = (i^3)*sum(EdMM2);

    % Calculate for Mindling and irb = 2 nad irs = 3
    
    ep = [ptype i 2 3];

    [AMM3, EdMM3 ] = deflectioncalc( ndofs,nelement, Edof, Ex, Ey,  DofTop, ...
        ep,eq,D,2,G,bc,fdofs,constraineddofs,tractiony,cdof,lx);

    deflectionMtotal3(counter,1) = (i^3)*sum(EdMM3);

    counter = counter +1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plot the moments
figure
subplot(3,1,1)
fill(Ex',Ey',Mxx');
colorbar
xlabel( colorbar,'Moment Distribution(N.m)');
title('Mxx')
subplot(3,1,2)
fill(Ex',Ey',Myy');
colorbar
xlabel( colorbar,'Moment Distribution(N.m)');
title('Myy')
subplot(3,1,3)
fill(Ex',Ey',Mxy');
colorbar
xlabel( colorbar,'Moment Distribution(N.m)');
title('Mxy')

% plot the stresses
figure
subplot(3,1,1)
fill(Ex',Ey',es(:,1)');
colorbar
xlabel( colorbar,'Stress Distribution(Pa)');
title('sigmaxx')
subplot(3,1,2)
fill(Ex',Ey',es(:,2)');
colorbar
xlabel( colorbar,'Stress Distribution(Pa)');
title('sigmayy')
subplot(3,1,3)
fill(Ex',Ey',es(:,3)');
colorbar
xlabel( colorbar,'Stress Distribution(Pa)');
title('sigmaxy')

% plot the highest compressive stress distribution
figure
fill(Ex',Ey',s_min')
colorbar
title('highest compressive stress')
xlabel( colorbar,'Stress Distribution(Pa)');
axis equal
axis([0 Lx 0 Ly])


% plot the Deflection versus thickness

thicknessvect = 0.001:0.01:0.5001;
figure


figure
plot(thicknessvect',deflectionMtotal1,'b-');
xlabel('Thickness(m)')
ylabel('Deflection(m)')
title('Thickness vs Deflection for Mindlin 2 1')
figure
plot(thicknessvect',deflectionMtotal2,'k-');
xlabel('Thickness(m)')
ylabel('Deflection(m)')
title('Thickness vs Deflection for Mindlin 2 2')
figure
plot(thicknessvect',deflectionMtotal3,'g-');
xlabel('Thickness(m)')
ylabel('Deflection(m)')
title('Thickness vs Deflection for Mindlin 2 3')
figure
plot(thicknessvect',deflectionKtotal,'r-');
xlabel('Thickness(m)')
ylabel('Deflection(m)')
title('Thickness vs Deflection for Kirchoff 2 1')

figure
hold on
plot(thicknessvect',deflectionMtotal1,'b-');
plot(thicknessvect',deflectionMtotal2,'k-');
plot(thicknessvect',deflectionMtotal3,'g-');
plot(thicknessvect',deflectionKtotal,'r-');

xlabel('Thickness(m)')
ylabel('Deflection(m)')
legend('Mindlin 2 1','Mindlin 2 2','Mindlin 2 3','Kirchoff 2 1')
title('Thickness vs Deflection')

hold off





