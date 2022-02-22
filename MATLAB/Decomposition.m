function [Dips,Strikes,Rakes,VOL,CLVD,DC,MDev]=Decomposition(M)

% This routine determines source parameters from moment tensor assuming the
% moment tensor M is in Carthesian coordinates defined by Aki and
% Richard, 1980, fig 4.20, page 114

Trace=trace(M)/3;
VOL=Trace/max(abs(eig(M)));
% After determination of the Volumetric part, substract it from the moment
% tensor and create Deviatoric moment tensor;
MDev=M-diag([Trace Trace Trace]);
% Determine eigen vectors and eigen values of the Deviatoric tensors,
[V,D]=eig(MDev);
[EDev,Order]=sort(diag(D));
% fault normal and displacement direction
n1=V(:,Order(1))+V(:,Order(3));
n1=n1/norm(n1,2);
u1=V(:,Order(1))-V(:,Order(3));
u1=u1/norm(u1,2);
% fault normal has negative z-component
if n1(3)>0
   n1=-n1;
end
% Correct for right sense of displacement
if (n1'*M*u1+u1'*M*n1)<0
    u1=-u1;
end
% Dips are angles between horizontal plane and fault plane
Dips(1)=acosd(abs(n1(3)));
% Strikes are angles from north (coordinate 1) and intersection of 
% the fault plane with horizontal plane 
Strikes(1)=180*atan2(-n1(1),n1(2))/pi;
if Strikes(1)<0
    Strikes(1)=360+Strikes(1);
end
% Rake is an angle between strike and slip
HStrike=[cosd(Strikes(1)) sind(Strikes(1)) 0];
Rakes(1)=acosd(dot(HStrike,u1));
if u1(3)>0
    Rakes(1)=-Rakes(1);
end

% Second Solution:
% fault normal and displacement direction
u1=V(:,Order(1))+V(:,Order(3));
u1=u1/norm(u1,2);
n1=V(:,Order(1))-V(:,Order(3));
n1=n1/norm(n1,2);
% fault normal has negative z-component
if n1(3)>0
    n1=-n1;
end
% Correct for right sense of displacement
if (n1'*M*u1+u1'*M*n1)<0
    u1=-u1;
end
% Dip is an angle between horizontal plane and fault plane
Dips(2)=acosd(abs(n1(3)));
% Strikes are angles from north (coordinate 1) and intersection of 
% the fault plane with horizontal plane 
Strikes(2)=180*atan2(-n1(1),n1(2))/pi;
if Strikes(2)<0
    Strikes(2)=360+Strikes(2);
end
% Rake is an angle between strike and slip
HStrike=[cosd(Strikes(2)) sind(Strikes(2)) 0];
Rakes(2)=acosd(dot(HStrike,u1));
if u1(3)>0
    Rakes(2)=-Rakes(2);
end

% Finally determine relative percentages of the DC and CLVD components
[Me,Ind]=min(abs(EDev));
Epsilon=EDev(Ind)/max(abs(EDev));
CLVD=-2*Epsilon*(1-abs(VOL));
DC=1-abs(CLVD)-abs(VOL);