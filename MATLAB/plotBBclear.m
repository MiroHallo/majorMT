function plotBB(MTs)
% Plot beach-ball into the current axes
% Input is the full Moment Tensor (1x6 vextor Aki & Richards notation)

hold on
Nc = length(MTs(:,1));

% -------------------------------------------------
% P and T quadrants
if Nc==1
    MT3 = [
        MTs(1,1) MTs(1,4) MTs(1,5)
        MTs(1,4) MTs(1,2) MTs(1,6)
        MTs(1,5) MTs(1,6) MTs(1,3)
        ];
    
    [x,y,z,c,~] = xyzc(MT3,'P',500);
    [X,Y,Z] = BProj(x,y,z);
    
    % Fill by colors
    maxc = max([abs(max(max(c))), abs(min(min(c)))]);
    colormap([1 1 1; 0 0 0])
    surf(X,Y,Z,c);
    shading interp;
    
    %Set the c scale to be symmetrical
    caxis([-maxc maxc]);
    set(gca,'CLim',[-maxc maxc])
    
    %Plot contours
    contour(X,Y,c,[0,0],'k');
else
    [x,y,z,c,~] = xyzc([-1 0 0;0 -1 0;0 0 -1],'P',100);
    [X,Y,Z] = BProj(x,y,z);
    
    % Fill by colors
    maxc = max([abs(max(max(c))), abs(min(min(c)))]);
    colormap([0.9 0.9 0.9; 0.9 0.9 0.9])
    surf(X,Y,Z,c);
    shading interp;
    
    %Set the c scale to be symmetrical
    caxis([-maxc maxc]);
    set(gca,'CLim',[-maxc maxc])
    
    %Plot contours
    contour(X,Y,c,[1,1],[0.9 0.9 0.9]);
end


% -------------------------------------------------
% Circle around the beach-ball
[x,y,z] = GetGreatCircle([1,0,0],[0,1,0]);
[X,Y,Z] = BProj(x,y,z);

plot3(X,Y,Z,'k');


% -------------------------------------------------
% DC nodal planes
for n=1:Nc
    MT3 = [
        MTs(n,1) MTs(n,4) MTs(n,5)
        MTs(n,4) MTs(n,2) MTs(n,6)
        MTs(n,5) MTs(n,6) MTs(n,3)
        ];
    
    [~,D] = eig(MT3);
    Lvect = [D(1,1), D(2,2), D(3,3)];
    
    if max(Lvect)<0 || min(Lvect>0)
        continue
    else
        [strike1,dip1,rake1,strike2,dip2,rake2] = MT2SDR(MT3);
        
        [vx1,vy1,vz1]=toa2xyz((strike1+90),(90-dip1),false,'deg');
        [vx2,vy2,vz2]=toa2xyz(strike1,90,false,'deg');
        [x1,y1,z1]=GetGreatCircle([vx1,vy1,vz1],[vx2,vy2,vz2]);
        
        [vx1,vy1,vz1]=toa2xyz((strike2+90),(90-dip2),false,'deg');
        [vx2,vy2,vz2]=toa2xyz(strike2,90,false,'deg');
        [x2,y2,z2]=GetGreatCircle([vx1,vy1,vz1],[vx2,vy2,vz2]);
        
        
        [X,Y,Z] = BProj(x1,y1,z1);
        plot3(X,Y,Z,'Color',[0.6 0.6 0.6],'LineWidth',1,'Marker','None')
        
        [X,Y,Z] = BProj(x2,y2,z2);
        plot3(X,Y,Z,'Color',[0.6 0.6 0.6],'LineWidth',1,'Marker','None')
    end
end


% -------------------------------------------------
% PT-axes
for n=1:Nc
    MT3 = [
        MTs(n,1) MTs(n,4) MTs(n,5)
        MTs(n,4) MTs(n,2) MTs(n,6)
        MTs(n,5) MTs(n,6) MTs(n,3)
        ];
    [T,N,P,E] = MT2TNP(MT3);
    
    if Nc==1
        [X,Y,Z] = BProj([T(1),-T(1)],[T(2),-T(2)],[T(3),-T(3)]);
        ph(1) = plot3(X,Y,Z,'s','MarkerEdgeColor','b','MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.5);
        %text(X,Y,Z+0.1,'  T','FontSize',10,'color','b');
        [X,Y,Z] = BProj([P(1),-P(1)],[P(2),-P(2)],[P(3),-P(3)]);
        ph(2) = plot3(X,Y,Z,'x','MarkerEdgeColor','r','MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.5);
        %text(X,Y,Z+0.1,'  P','FontSize',10,'color','b');
    else
        [X,Y,Z] = BProj([T(1),-T(1)],[T(2),-T(2)],[T(3),-T(3)]);
        ph(1) = plot3(X,Y,Z,'s','MarkerEdgeColor','b','MarkerSize',10,'LineWidth',1.5);
        [X,Y,Z] = BProj([P(1),-P(1)],[P(2),-P(2)],[P(3),-P(3)]);
        ph(2) = plot3(X,Y,Z,'x','MarkerEdgeColor','r','MarkerFaceColor','w','MarkerSize',10,'LineWidth',1.5);
    end
end


% -------------------------------------------------
% Final formating
axis normal;
axis equal;
axis off;
grid off;
set(gca,'view',[-90 90])
set(gca,'Ydir','rev');
%legend(ph,'T-axis','P-axis','Location','southoutside')

hold off

end


% -------------------------------------------------
% -------------------------------------------------
% -------------------------------------------------
% help fuction for projection (lower)
function [X,Y,Z]=BProj(x,y,z)
X=x.*sqrt(2./(1+z));
Y=y.*sqrt(2./(1+z));
X(z<0)=NaN;
Y(z<0)=NaN;
Z=zeros(size(X));
X(isinf(X))=NaN;
Y(isinf(Y))=NaN;
end


