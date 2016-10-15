function [AMM, Ed ] = deflectioncalc( ndofs,nelement, Edof, Ex, Ey,  DofTop, ...
    ep,eq,D,plate,G,bc,fdofs,constraineddofs,tractiony,cdof,Lx )
% function [AMM, Ed ] = deflectioncalc( ndofs,nelement, Edof, Ex, Ey,  DofTop, 
% ep,eq,D,plate,G,bc,fdofs,constraineddofs,tractiony,cdof,Lx )
% 
% -------------------------------------------------------------------------
% Purpose : This function will calculate the Deflection for a give thickness
% and stress for the centermost degrees of freedom (cdof)
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
% plate                     Kirchoff = 1 Mindlin = 2
%
% ------------------------------------------------------------------------
% Output
%
% Ed = extract(Edof,AMM); 
% AMM = ndofs x 1           Displacements


[AM,~,~]  = displacementcalc( ndofs,nelement,Edof, Ex, Ey,  DofTop, ...
    ep,eq,D,plate,G,bc,fdofs,constraineddofs,tractiony*ep(2),Lx);

AMM = zeros(ndofs,1);
AMM(cdof) = AM(cdof);

Ed = extract(Edof,AMM);

Ed = Ed(:,3:5:end);

% Since I only calculated the displacements for cdof
poszeros  = find(Ed == 0);

Ed(poszeros) = [];

end

