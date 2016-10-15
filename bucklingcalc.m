function [ X,L,posvec,freedof,realambdaoriginal ] = bucklingcalc( Ex,Ey,ep,D,es,G...
                                                    ,plate2nd,eq,bc,ndofs,Edof,nelement )
% function[ X,L,posvec,freedof,realambdaoriginal ] = bucklingcalc( Ex,Ey,ep,D,es,G...
%          ,plate2nd,eq,bc,ndofs,Edof,nelement )
% -------------------------------------------------------------------------
% Purpose : This function will calculate the eingenvalues L, the realeingenvalues
% and the matrix X representes the eigenmodes. It will give the position of 
% the real eingenvalues , ordered from the lowest to the biggest eigenvalue, 
% in the posvec vector.
%
% -------------------------------------------------------------------------
% Input
%
% ndofs                     Total number of degrees of freedom
% nelement                  Total number of elements
% Ex                        Element nodal x-coordinates, where the nodes are
%                           counted anti-clockwise starting from the bottom left
%                           node of each element.
% Ey                        Element nodal y-coordinates, where the nodes are
%                           counted anti-clockwise starting from the bottom left
%                           node of each element.
% bc                        Boundary Conditions
% G    
% es = [sigxx sigyy sigxy] stresses
% eq = [qx qy qz]               
% ep = [ptype t irb irs]    ptype: analysis type,    t: thickness                                 
% D                         Hook matrix   
% plqte2nd = 1              1 for Kirchoff
%
% ------------------------------------------------------------------------
% Output
%
% X                  eingenvectors 
% L                  enigenvalues
% posvec             vector with the positions of the real eigenvectors 
%                    ordered from the smallest to the biggest
% realambdaoriginal  real eigenvectors
% freedof            free degrees of freedom 




KM = zeros(ndofs,ndofs);
KG = zeros(ndofs,ndofs);
FB = zeros(ndofs, 1);
% Compute KM and KG for each element and assemble them

for gamma = 1:nelement
    
    elemdofs = Edof(gamma,2:end);
    
    [KMe,fe,KGe] = shell2rg(Ex(gamma,:),Ey(gamma,:),ep,D,es(gamma,:),G,plate2nd,eq);
    
    KM(elemdofs,elemdofs) = KM(elemdofs,elemdofs) + KMe;
    
    KG(elemdofs,elemdofs)= KG(elemdofs,elemdofs) +KGe;
    
    FB(elemdofs,1) = FB(elemdofs,1) + fe;
    
end


% Solving the eingenvalue problem

bc1 = bc(:,1);

freedof = setdiff(1:ndofs,bc1);

Kred = sparse(KM(freedof,freedof));
Gred = sparse(KG(freedof,freedof));

Knew = Kred\Gred;

[X,L] = eigs(Knew);

realambda = -1./diag(L)

realambdaoriginal = zeros(1,length(realambda));
posvec = zeros(length(realambda),1);


% Recovering the position and ordering it from the minimum to the biggest
% lambda

for u = 1:length(realambda)
    
    lambdamin = min(abs(realambda));
    posmin = find(realambda == lambdamin);
    realambdaoriginal(1,u) = realambda(posmin);
    realambda(posmin) = [];
    posvec(u,1) = posmin + (u-1);
end

end

