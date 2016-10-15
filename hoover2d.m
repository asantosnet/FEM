function [Ex,Ey,Edof,Boundarygeom,Boundarydof]=hoover2d(nthick,nheight)
% [Ex,Ey,Edof,Boundarygeom,Boundarydof]=hoover2D(nthick,tneight)
%
% Purpose:  Generates the mesh of a 2D cross section of the Hoover dam
%           consisting of bilinear quadrilateral elements.
%
% Input:    nthick       - number of elements in the thickness direction
%           nheight      - numbet of elements in the hight direction
% Output:   Ex           - Element x-coordinates
%           Ey           - Element y-coordinates
%           Edof         - Element degrees of freedom
%           Boundarygeom - Boundary segment geometric data
%               Boundarygeom=[x_1^1 y_1^1 x_2^1 y_2^1 n_x^1 n_y^1;
%                               .     .     .     .     .     .  ;
%                             x_1^n y_1^n x_2^n y_2^n n_x^1 n_y^n]
%               Data row f,
%               x_1^f    - x-coordinate, node 1, edge f
%               y_1^f    - y-coordinate, node 1, edge f
%               x_2^f    - x-coordinate, node 2, edge f
%               y_2^f    - y-coordinate, node 2, edge f
%               n_x^f    - x-component of (outwars) normal vector, edge f
%               n_y^f    - y-component of (outwars) normal vector, edge f
%           Boundarydof  - Boundary segment degrees of freedom
%               Boundarydof=[dof_x1^1 dof_y1^1 dof_x2^1 dof_y2^1;
%                              .     .     .     .     .     .  ;
%                            dof_x1^n dof_y1^n dof_x2^n dof_y2^n]
%               Data row f,
%               dof_x1^n - x-degree of freedom, node 1, edge f
%               dof_y1^n - y-degree of freedom, node 1, edge f
%               dof_x2^n - x-degree of freedom, node 2, edge f
%               dof_y2^n - y-degree of freedom, node 2, edge f

close all

xi=[];
eta=[];

for i=1:(nheight+1)
    xi=[xi linspace(0,1,nthick+1)];
    eta=[eta ones(1,nthick+1).*(i-1)/(nheight)];
    Enod(nthick*(i-1)+(1:nthick),1:4)=[((nthick+1)*(i-1)+(1:nthick))' ((nthick+1)*(i-1)+(2:(nthick+1)))' ((nthick+1)*i+(2:(nthick+1)))' ((nthick+1)*i+(1:nthick))'];
end

Bnod(1:nthick,:)=[(1:(nthick))' (2:(nthick)+1)'];
Bnod(nthick+(1:nheight),:)=(nthick+1).*[(1:(nheight))' (2:(nheight+1))'];
s=1:(nthick+1);
Bnod(nthick+nheight+(1:nthick),:)=(nthick+1)*(nheight+1)+1-[s(1:(end-1))' s(2:end)'];
s=1:(nheight+1);
Bnod(2*nthick+nheight+(1:nheight),:)=(nthick+1)*(nheight)+1+(nthick+1).*(1-[s(1:(end-1))' s(2:end)']);

nel=nthick*nheight;
nnod=(nthick+1).*(nheight+1);
Enod=Enod(1:nel,:);

%Trans-finite interpolation
y=eta.*724;

x=(1-xi).*leftside(y)+xi.*rightside(y);

x=x.*0.3048;
y=y.*0.3048;

Ex=x(Enod)
Ey=y(Enod)

eldraw2(Ex,Ey,[1 1 0]);

hold on
plot([0 0],[0 250],'-.b')
xlabel('x=r [m]')
ylabel('y [m]')
plot([rightside(200/.3048).*.3048 350],[200 200])

%axis([-10 350 0 250])
text(300,210,'y=200 m')

Boundarygeom(:,1:4)=[x(Bnod(:,1))' y(Bnod(:,1))' x(Bnod(:,2))' y(Bnod(:,2))'];
Boundarygeom(:,5:6)=[y(Bnod(:,2))'-y(Bnod(:,1))' -(x(Bnod(:,2))'-x(Bnod(:,1))')]./(sqrt((x(Bnod(:,2))'-x(Bnod(:,1))').^2+(y(Bnod(:,2))'-y(Bnod(:,1))').^2)*ones(1,2));

Dof=[(1:2:(2*nnod-1))' (2:2:(2*nnod))'];

Edof(:,1)=(1:nel)';
Edof(:,2:2:8)=[Dof(Enod(:,1),1) Dof(Enod(:,2),1) Dof(Enod(:,3),1) Dof(Enod(:,4),1)];
Edof(:,3:2:9)=[Dof(Enod(:,1),2) Dof(Enod(:,2),2) Dof(Enod(:,3),2) Dof(Enod(:,4),2)];

Boundarydof=[Dof(Bnod(:,1),1) Dof(Bnod(:,1),2) Dof(Bnod(:,2),1) Dof(Bnod(:,2),2)];

scl=10;
for i=1:length(Boundarygeom(:,1))
    plot((Boundarygeom(i,1)+Boundarygeom(i,3))/2+[0 scl].*Boundarygeom(i,5),(Boundarygeom(i,2)+Boundarygeom(i,4))/2+[0 scl].*Boundarygeom(i,6))
end

function r=leftside(z)

k1=412/482;
k2=90/242;

r=218+k1.*z+((k2-k1).*(z-482)).*signp(z-482);

function r=rightside(z)

A=(878-765)/(724^2);

r=A.*(724-z).^2+765;

function y=signp(x)

y=(sign(x)+1)/2;