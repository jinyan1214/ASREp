function [N,dNdx,J] = shape3d (nat,xyz)
% SHAPE3D shape functions for 8 node solid (brick) element
% [N, dNdx, J] = SHAPE3D (NAT,XYZ) shape functions for 8 node solid (brick) element
%%
% Input Parameters
% ----------------
% nat  = [ xi eta zeta ] natural coordinates in element
% xyz  = nodal coordinates for element (row i for node i)
% ----------------
% Return Variables
% ----------------
% N    = shape functions in natural coordinates
% dNdx = dNdx(i,j) = derivative of shape function j with respect to geometric coordinate x_i
% J    = Jacobian of transformation from geometric to natural coordinates
% ------------------------------------------------------

%  =========================================================================================
%  FEDEASLab - Release 3.0, July 2008
%  Matlab Finite Elements for Design, Evaluation and Analysis of Structures
%  Professor Filip C. Filippou (filippou@ce.berkeley.edu)
%  Department of Civil and Environmental Engineering, UC Berkeley
%  Copyright(c) 1998-2008. The Regents of the University of California. All Rights Reserved.
%  =========================================================================================

%%
% Nodal coordinates in natural coordinate system
natnode = [ -1 +1 +1 -1 -1 +1 +1 -1;
            -1 -1 +1 +1 -1 -1 +1 +1;
            -1 -1 -1 -1 +1 +1 +1 +1];

% terms in shape functions
xi   = 0.5 * ( 1 + natnode(1,:) * nat(1) );
eta  = 0.5 * ( 1 + natnode(2,:) * nat(2) );
zeta = 0.5 * ( 1 + natnode(3,:) * nat(3) );

% shape functions
N = xi .* eta .* zeta;

% derivatives of shape functions with respect to natural coordinates
dN(1,:) = 0.5 * natnode(1,:) .* eta .* zeta;   % dN/d-xi
dN(2,:) = 0.5 * natnode(2,:) .* xi  .* zeta;   % dN/d-eta
dN(3,:) = 0.5 * natnode(3,:) .* xi  .* eta;    % dN/d-eta

% compute Jacobian of coordinate transformation
J    = dN * xyz;
% compute derivatives of shape functions with respect to geometric coordinates
dNdx = J \ dN;