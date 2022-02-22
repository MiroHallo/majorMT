function [Strike1,Dip1,Rake1,Strike2,Dip2,Rake2] = MT2SDR(MT3)
% Calculate the strike dip and rake of the nodal planes from the MT
% Using MT2TNP, TP2FP and FP2SDR
% Input is the full MT (3x3 matrix Aki & Richards notation)
% Coordinates are x=North, y=East, z=Down
%
% Converts a moment tensor to the T, N, P axes
[T,N,P,E] = MT2TNP(MT3);

% Calculate the normals to the nodal planes from the T and P axes
[N1,N2] = TP2FP(T,P);

% Calculate strike, dip and rake for a normal and slip vector
[Strike1,Dip1,Rake1] = FP2SDR(N1,N2);
[Strike2,Dip2,Rake2] = FP2SDR(N2,N1);
end
