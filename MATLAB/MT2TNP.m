function [T,N,P,E] = MT2TNP(MT3)
% Converts a moment tensor to the T, N, P axes
% Input is the full Moment Tensor (3x3 matrix Aki & Richards notation)
% Returns the T, N, P eigen vectors and the eigenvalues in E as (v1,v2,v3)
%
% Remove ISO
% Trace = trace(MT3)/3;
% ISO = Trace/max(abs(eig(MT3)));
% MT3 = MT3-diag([Trace Trace Trace]);

% Eigenvalues and eigenvectors
[V,D] = eig(MT3);
Lvect = [D(1,1), D(2,2), D(3,3)];
[L,Li] = sort(Lvect,'descend');
Et = D(Li(1),Li(1)); % Lambla 1
En = D(Li(2),Li(2)); % Lambla 2
Ep = D(Li(3),Li(3)); % Lambla 3
T = V(:,Li(1)); % Eigenvector 1
N = V(:,Li(2)); % Eigenvector 2
P = V(:,Li(3)); % Eigenvector 3
E = [Et En Ep];

% Correct P to be pointing downwards -> normals oriented correctly
if P(3)<0
    P = -P;
end
end