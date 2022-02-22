%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Decompose non-DC MT into Major and Minor DC MTs with preserved dominant
% eigenvector (P- or T-axis)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Method by Hallo et al. (2017, 2019)
% Hallo,M., Asano,K., Gallovic,F. (2017): Bayesian inference and interpretation
%      of centroid moment tensors of the 2016 Kumamoto earthquake sequence, 
%      Kyushu, Japan, Earth, Planets and Space, 69:134.
% Hallo,M., Oprsal,I., Asano,K., Gallovic,F. (2019): Seismotectonics of the 2018
%      Northern Osaka M6.1 earthquake and its aftershocks: joint movements
%      on strike-slip and reverse faults in inland Japan, Earth,
%      Planets and Space, 71:34.
%
% Code author: Miroslav Hallo
% Charles University in Prague, Faculty of Mathematics and Physics
% Web: http://geo.mff.cuni.cz/~hallo/
% E-mail: hallo@karel.troja.mff.cuni.cz
% Revision 2/2017: The first version of the function.
% Revision 12/2018: Enhanced version.
%
% Copyright (C) 2017,2018  Miroslav Hallo
%
% This program is published under the GNU General Public License (GNU GPL).
%
% This program is free software: you can modify it and/or redistribute it
% or any derivative version under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
% and don't remove their names from the code.
%
% You should have received copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;
addpath([pwd,'/MATLAB']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT
% Moment tensor in Harvard notation
%MT = [-1.66  5.06  -3.40  0.375  0.94  0.824]; % 2016 Kumamoto M7.2
%MT = [-0.42  1.56  -1.14  -0.33  -0.53  0.73]; % 2016 Kumamoto M6.5
%MT = [-4.26  4.55  -0.29  -3.16  -1.46  1.69]; % 2016 Kumamoto M5.4
MT = [1.10  1.53  -2.65  0.26  0.08  -0.72]; % 2018 Osaka 6.1

% END INPUT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[MTa] = hrv2ar(MT);

% Figure of full input MT
figure('color','w')
plotBB(MTa);
title('Input full MT')
axis equal
axis off

% From Harvard to Aki & Richards notation
MT3 = [
MTa(1) MTa(4) MTa(5)
MTa(4) MTa(2) MTa(6)
MTa(5) MTa(6) MTa(3)
];

% Remove ISO
Trace = trace(MT3)/3;
ISO = Trace/max(abs(eig(MT3)));
MT3 = MT3-diag([Trace Trace Trace]);

% Figure of input deviatoric MT
figure('color','w')
plotBB([MT3(1,1),MT3(2,2),MT3(3,3),MT3(1,2),MT3(1,3),MT3(2,3)]);
title('Input deviatoric MT')
axis equal
axis off

% Eigenvalues and eigenvectors
[V,D] = eig(MT3);

% Sort eigenvalues and eigenvectors
Lvect = [D(1,1), D(2,2), D(3,3)];
display('I--------------------------------')
display('Non-DC MT eigenvalues:')
display(num2str(Lvect))
display(' ')
[L,Li] = sort(abs(Lvect),'descend');
La1 = D(Li(1),Li(1)); % Lambla 1
La2 = D(Li(2),Li(2)); % Lambla 2
La3 = D(Li(3),Li(3)); % Lambla 3
v1 = V(:,Li(1)); % Eigenvector 1
v2 = V(:,Li(2)); % Eigenvector 2
v3 = V(:,Li(3)); % Eigenvector 3

% Major and minor DC MTs
MTmaj = La2*(- v1*v1' + v2*v2'); % Major DC 
MTmin = La3*(- v1*v1' + v3*v3'); % Minor DC

% Dispay decomposed MT in Harvard notation
display('II-------------------------------')
display('Input non-DC MT (Harvard):')
[mt_tmp(1),mt_tmp(2),mt_tmp(3),mt_tmp(4),mt_tmp(5),mt_tmp(6)] = ar2hrv(MT3);
display(num2str(mt_tmp))
display(' ')
display('Major DC MT (Harvard):')
[mt_tmp(1),mt_tmp(2),mt_tmp(3),mt_tmp(4),mt_tmp(5),mt_tmp(6)] = ar2hrv(MTmaj);
display(num2str(mt_tmp))
display(' ')
display('Minor DC MT (Harvard):')
[mt_tmp(1),mt_tmp(2),mt_tmp(3),mt_tmp(4),mt_tmp(5),mt_tmp(6)] = ar2hrv(MTmin);
display(num2str(mt_tmp))
display(' ')

% Decompose Major and Minor MT
[Dips1,Strikes1,Rakes1,VOL1,CLVD1,DC1,MDev1] = Decomposition(MTmaj);
[Dips2,Strikes2,Rakes2,VOL2,CLVD2,DC2,MDev2] = Decomposition(MTmin);

display('III-----------------------------')
display('Major DC (Strike/Dip/Rake):')
display(num2str([Strikes1(1), Dips1(1), Rakes1(1)],'%6.1f '))
display(num2str([Strikes1(2), Dips1(2), Rakes1(2)],'%6.1f '))
display(' ')
display('Minor DC (Strike/Dip/Rake):')
display(num2str([Strikes2(1), Dips2(1), Rakes2(1)],'%6.1f '))
display(num2str([Strikes2(2), Dips2(2), Rakes2(2)],'%6.1f '))
display(' ')
display('DC percentage of Major and Minor DC MTs:')
display([num2str(DC1*100),'% vs ',num2str(DC2*100),'%'] )
display(' ')

% Ratio of scalar moments of Major and Minor DC MTs
Cdc1 = max(abs(eig(MTmaj)));
Cdc2 = max(abs(eig(MTmin)));
N = Cdc1+Cdc2;
Cdc1 = Cdc1 / N;
Cdc2 = Cdc2 / N;
display('Scalar seismic moment ratio of Major and Minor DC MTs:')
display([num2str(Cdc1*100,'%5.1f'),'% vs ',num2str(Cdc2*100,'%5.1f'),'%'] )

% Plot Major and Minor MT
figure('color','w')
subplot(1,3,1)
bb([Strikes1(1) Dips1(1) Rakes1(1)], 0.5, 0.5, Cdc1/2, 0, 'k')
title('Major DC')
axis equal
set(gca,'XLim',[0 1])
set(gca,'YLim',[0 1])
axis off

subplot(1,3,2)
bb([Strikes2(1) Dips2(1) Rakes2(1)], 0.5, 0.5, Cdc2/2, 0, 'k')
title('Minor DC')
axis equal
set(gca,'XLim',[0 1])
set(gca,'YLim',[0 1])
axis off

subplot(1,3,3)
plotBBclear([MT3(1,1),MT3(2,2),MT3(3,3),MT3(1,2),MT3(1,3),MT3(2,3)]);
title('Deviatoric non-DC')
axis equal

axis off





return
