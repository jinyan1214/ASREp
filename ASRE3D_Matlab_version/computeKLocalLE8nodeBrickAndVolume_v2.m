function [klocal, vol] = computeKLocalLE8nodeBrickAndVolume_v2(xyzLocal, E, v, G)


% 2x2x2 gauss integration points
nat  = [ -1 -1 -1; +1 -1 -1; +1 +1 -1; -1 +1 -1;
               -1 -1 +1; +1 -1 +1; +1 +1 +1; -1 +1 +1] * (1/sqrt(3));
wIP = ones(1,8);
nIP = length(wIP);
natnode = [ -1 +1 +1 -1 -1 +1 +1 -1;
            -1 -1 +1 +1 -1 -1 +1 +1;
            -1 -1 -1 -1 +1 +1 +1 +1];
klocal = zeros(24, 24);
% lamda = E*v/(1+v)/(1-2*v);
% mu = E/2/(1+v);
% D = [lamda + 2*mu, lamda, lamda, 0, 0, 0;
%      lamda, lamda + 2*mu, lamda, 0, 0, 0;
%      lamda, lamda, lamda + 2*mu, 0, 0, 0;
%      0, 0, 0, mu, 0, 0;
%      0 ,0, 0, 0, mu, 0;
%      0, 0, 0, 0, 0, mu];
D = [1, v/(1-v), v/(1-v), 0, 0, 0;
        v/(1-v), 1, v/(1-v), 0, 0, 0;
        v/(1-v), v/(1-v), 1, 0, 0, 0;
        0, 0, 0, (1-2*v)/2/(1-v), 0, 0;
        0, 0, 0, 0, (1-2*v)/2/(1-v), 0;
        0, 0, 0, 0, 0, (1-2*v)/2/(1-v)] * E * (1-v)/(1+v)/(1-2*v);
% G = E/12.5;

D(4,4) = G;
D(5,5) = G;
D(6,6) = G;
vol = 0;
for i = 1:nIP
    % terms in shape functions evaluated at gaussian points
    xi   = 0.5 * ( 1 + natnode(1,:) * nat(i, 1) );
    eta  = 0.5 * ( 1 + natnode(2,:) * nat(i, 2) );
    zeta = 0.5 * ( 1 + natnode(3,:) * nat(i, 3) );
    
    % derivatives of shape functions with respect to natural coordinates
    % evaluated at gaussian points
    dN(1,:) = 0.5 * natnode(1,:) .* eta .* zeta;   % dN/d-xi
    dN(2,:) = 0.5 * natnode(2,:) .* xi  .* zeta;   % dN/d-eta
    dN(3,:) = 0.5 * natnode(3,:) .* xi  .* eta;    % dN/d-zeta
    
    Jmat = dN * xyzLocal; %determinate (CE222 note page 248)
    
    dNdx = (Jmat) \  dN; % Jmat^-1 * dN
    B = zeros(6, 24);
%     B([1 4 5],1:3:24) = dNdx;
%     B([4 2 6],2:3:24) = dNdx;
%     B([5 6 3],3:3:24) = dNdx;
    B([1 4 6],1:3:24) = dNdx;
    B([4 2 5],2:3:24) = dNdx;
    B([6 5 3],3:3:24) = dNdx;
    klocal = klocal + B' * D * B * det(Jmat)*wIP(i);
    vol = vol + det(Jmat)*wIP(i);
    
end
% make sure klocal symmetric
klocal = 0.5*(klocal + klocal'); 

% klocal = ones(24, 24);

return

