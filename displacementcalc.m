function [ A, Ed ,F] = displacementcalc(ndofs,nelement, Edof, Ex, Ey,  DofTop, ...
                                      ep,eq,D,plate,G,bc,fdofs,constraineddofs,tractiony,Lx)
% function [ A, Ed ,F] = displacementcalc(ndofs,nelement, Edof, Ex, Ey,  DofTop, ...
% ep,eq,D,plate,G,bc,fdofs,constraineddofs,tractiony,Lx)
% 
% -------------------------------------------------------------------------
% Purpose : This function will calculate the Displacement for a give
% thickness and stress. Mindlin and Kirchoff methods can be used
%
% -------------------------------------------------------------------------
% Input
%
% ndofs                     Total number of degrees of freedom
% nelement                  Total number of elements
% Edof                      Topology matrix in terms of dofs, size = (NoElem x 21)
%                           In each row in Edof, the dofs are arranged 
%                           anti-clockwise starting from the bottom left node
%                           of each element (5 dofs/node). For each node, first
%                           the three translational dofs are counted (ux,uy,uz) 
%                           and then the two rotational (theta_x,theta_y)
% Ex                        Element nodal x-coordinates, where the nodes are
%                           counted anti-clockwise starting from the bottom left
%                           node of each element.
% Ey                        Element nodal y-coordinates, where the nodes are
%                           counted anti-clockwise starting from the bottom left
%                           node of each element.
% DofTop                    Translational dofs (ux,uy,uz) at the top boundary,
%                           size = (NoNodesRightBoundary x 3)
% fdofs                     Free degrees of freedom
% bc                        Boundary Conditions
% G               
% constraineddofs           Contrained degrees of freedom
% eq = [qx qy qz]               
% ep = [ptype t irb irs]    ptype: analysis type,    t: thickness                                 
% D                         Hook matrix
% tractiony                 Traction along y direction
% cdof                      Dofs at the centermost point
% Lx                        Length along x of each element   
%
% ------------------------------------------------------------------------
% Output
%
% Ed = extract(Edof,AMM); 
% A = ndofs x 1           Displacements
% F

A = zeros(ndofs,1);

K = zeros(ndofs,ndofs);
F = zeros(ndofs,1);

for n = 1:nelement
    
    elemdofs = Edof(n,2:end);
    [Ke,Fe ] = shell2re(Ex(n,:),Ey(n,:),ep,D,G,eq,plate);
    
    K(elemdofs,elemdofs) = K(elemdofs,elemdofs) + Ke;
    F(elemdofs) = F(elemdofs) + Fe;
    
end


% tractionydof = (tractiony)/length(DofTop);

F(DofTop(5:3:(end-3),1)) = F(DofTop(5:3:(end-3),1))  + ((tractiony*(Lx)));
F(DofTop(2:(end-3):end,1)) = F(DofTop(2:(end-3):end,1)) + ((tractiony*Lx)/2);


A(fdofs) = (K(fdofs,fdofs))\(F(fdofs) - K(fdofs,constraineddofs)*bc(:,2));

Ed = extract(Edof,A);

end

