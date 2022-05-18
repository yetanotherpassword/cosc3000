clear all
close all


%fname="~/3pm_20jul-elev5.csv" 
fname="~/3pm_31dec-elev5.csv" 
opts = detectImportOptions(fname,'Delimiter',',','PartialFieldRule','keep','VariableNamingRule','preserve');
%opts.SelectedVariableNames = opts.SelectedVariableNames([9, 27, 29, 31, 32, 34, 38, 39, 44, 79]);


T = readtable(fname,opts);

%polarbubblechart(T,"RFS_Azimuth","RFS_Range",[],"RangeError_BT_2_and_10",'MarkerFaceAlpha',0.10);

%T.Properties.VariableNames=varNames;
% clear all
%   r=0:1:60;
%    theta=0:pi/30:2*pi;
%    for i=1:61
%               r(i)=i-1;
%               if (r(i)<3)
%                   sigmar(i)=-6;
%               else
%                   sigmar(i)=((1-((7^2)/(r(i)^2)))*-14);   %multiplier is a predefined constant
%               end
%    end
%     [x,y]=pol2cart(theta, r);
%     [X,Y]=meshgrid(x,y);
%     sigmar=repmat(sigmar,61,1);
%    surf(X,Y,sigmar)
qmin=min(T.RangeError_BT_2_and_10)+0.0001
siz=size(T,1);
a=[];
b=[];
c=[];
d=[];
id=[];
figure
for i=1:1:siz
    fprintf("%d\n",i)
    if isnan(T.RFS_Azimuth(i))
        c=c+qmin
        %m=max(c)
        %c1=round(c/m*6);
        %c2=ones(size(c1,1))*2;
        %c3=c2.^c1;
        id=[id; T.MSID(i+1)]
        c3=c;
        polarbubblechart(a,b,c3,'MarkerFaceAlpha',0.05);
  
        hold on
        a=[];
        b=[];
        c=[];
        d=[];
    else
        a=[ a; T.RFS_Azimuth(i)*pi/180 ];
        b=[ b; T.RFS_Range(i)];
        c=[ c; T.RangeError_BT_2_and_10(i)];
        d=[d; T.Elevation(i)]
    end
end

bubblelegend('Range Error (m)','Location','eastoutside')
legend(id)
title("Polar Plot of Individual Tracks with Range Error > {\sigma}")
subtitle("Mt Macedon, Melbourne 31/12/2021")
rmax = 256;
%hax = polaraxes('RLim', [0 rmax]);
text(0, -rmax/2.5, 'Range from Radar (NM)', 'horiz', 'center', 'vert', 'top', 'rotation', 0);
text(3*pi/4, rmax/1.5, 'Bearing from Radar (0^o=NORTH)', 'horiz', 'center', 'rotation', 45);
function polarscatter3(theta,phi,rho, varargin )
    maxr=1;
    if nargin>2, maxr = max(ceil(max(rho)),1); end

    %% Plot axis
    clrh = 0;
    if ~ishold, clf, clrh = 1; end 
    hold on;
    for r = maxr/5:maxr/5:maxr
        [tax_x,tks_y] = pol2cart(linspace(0,2*pi,101),r);
        tax_z = zeros(length(tax_x),1); 
        pax_x = tax_z; pax_y = tax_x; pax_z = tks_y;
                
        if r==maxr
            plot3(pax_x,pax_y,pax_z,'-','Color','#606060','LineWidth',1)
            patch(pax_x,pax_y,pax_z,'w','FaceAlpha',0.6)
            plot3(tax_x,tks_y,tax_z,'-','Color','#606060','LineWidth',1)
            patch(tax_x,tks_y,tax_z,'w','FaceAlpha',0.6)
        else
            plot3(pax_x,pax_y,pax_z,'-','Color','#AAAAAA','LineWidth',0.1)
            plot3(tax_x,tks_y,tax_z,'-','Color','#AAAAAA','LineWidth',0.1) 
        end
    end
    plot3(-1:1,[0,0,0],[0,0,0],'-','Color','#000000','LineWidth',1)
    plot3([0,0,0],-1:1,[0,0,0],'-','Color','#000000','LineWidth',1)
    plot3([0,0,0],[0,0,0],-1:1,'-','Color','#000000','LineWidth',1)
    
    lim = maxr+0.1;
    xlim([-lim lim]), ylim([-lim lim]), zlim([-lim lim])

    
    set(gca,'xtick',[]),set(gca,'ytick',[]),set(gca,'ztick',[]);
    axis off;
    view([135,45])

    %% Plot ticks
        r = 1;
     for d = 0:pi/10:2*pi
        [tks_x,tks_y] = pol2cart([d d],[r r+0.05]); 
        cero = [0 0];
        plot3(tks_x,tks_y,cero,'-','Color','#606060','LineWidth',1)
        plot3(cero,tks_x,tks_y,'-','Color','#606060','LineWidth',1)
     end
     
     %%Plot scater of points
     if nargin <2
         return;
     end
     [x,y,z] = sph2cart(theta,phi,rho);
     if nargin<4
         scatter3(x,y,z)
     else
         scatter3(x,y,z,varargin{:})
     end
     
     if clrh, hold off; end
end

