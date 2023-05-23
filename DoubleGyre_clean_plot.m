% plot the Model and GPR predicted doubel gyre velocity fields
% Just in case you are not familar with matlab code, I simlified the code
% into 2 case, changing the 'Casenumber' and see the difference between the
% resluts.

% Basic settings, stream function 
% phi = sin(x)*sin(y)+epsilon*sin(x)-omega*t*sin(2*y);
% epsilon=0.1; omega=2*pi/10; according to previous paper
clear;clc;

Casenumber=2;
% Casenumber=1 for a longer time interval ts=0.5
% Casenumber=2 for a shorter time interval ts=0.2


switch Casenumber
    case 1
        % Below is dataset1 for a longer time interval 0.5 load afterwards
%         load('DoubleGyre_50drifter_40timestep_DD_OI.mat');
        ts=0.5
    case 2
        % Below is dataset2 for a shorter time interval 0.2 load afterwards
%         load('DoubleGyre_50drifter_100timestep_DD_OI.mat');
        ts=0.2     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Two Gyre Model for velocity fields %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon=0.1;
omega=2*pi/10;
dxy=0.05;
x=[0:dxy:6+29*0.01];y=[0:dxy:3+14*0.01];
I=length(x);J=length(y);
[MX,MY]=meshgrid(x,y);
% t=[0:ts:20-0.1];
t=[0:0.2:20];
SkipX=0.2/dxy; % for the space of arrow
SkipY=0.2/dxy;

for ti=1:length(t)
    for i=1:I
        for j=1:J
            % for stream function should be u=phi / y; v=- phi/ x; might be
            % some sign problem
            phi(j,i)=sin(x(i))*sin(y(j))+epsilon*sin(x(i)-omega*t(ti))*sin(2*y(j));
            u(j,i)=-(sin(x(i))*cos(y(j))+epsilon*sin(x(i)-omega*t(ti))*2*cos(2*y(j)));
            v(j,i)=cos(x(i)).*sin(y(j))+epsilon*cos(x(i)-omega*t(ti)).*sin(2*y(j));
            
        end
    end
    U(ti)={u};V(ti)={v};
    figure(1)
    clf(1)
%     l=streamslice(MX(1:SkipX:end,1:SkipY:end),MY(1:SkipX:end,1:SkipY:end),u(1:SkipX:end,1:SkipY:end),v(1:SkipX:end,1:SkipY:end));
%     set(l,'color','b')
    
    box on;
    contour(x,y,phi);
    hold on;
    quiver(MX(1:SkipX:end,1:SkipY:end),MY(1:SkipX:end,1:SkipY:end),u(1:SkipX:end,1:SkipY:end),v(1:SkipX:end,1:SkipY:end));
    axis([0,6.29,0,3.14])
    saveas(figure(1),['DoubleGyre_TimeStep',num2str(ti),'.png'])
    hold off;
%     pause(0.3)
    clear phi u v
end
axis([0,6.29,0,3.14])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Advecting Lagrangian Particles %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set drifter locations
n=50;
T_limit=[0 0];
X_limit=[1 5];Y_limit=[0.5 2.5];
vel_u=U{1}';
vel_v=V{1}';

% Drifter Advection
ID=[1:n]';
% fix random
rng('default')
s=rng;
% location and time dispersion
for i=1:n
    T(i,:)=floor(T_limit(1)+(T_limit(2)-T_limit(1))*rand);
    X(i,:)=X_limit(1)+(X_limit(2)-X_limit(1))*rand;
    Y(i,:)=Y_limit(1)+(Y_limit(2)-Y_limit(1))*rand;

% Velocity interpolation
    % find the nearest 16 points
    Xrange=[floor(X(i,:)/dxy)-1:ceil(X(i,:)/dxy)+1];
    Yrange=[floor(Y(i,:)/dxy)-1:ceil(Y(i,:)/dxy)+1];
    Xq=X(i,:);Yq=Y(i,:);
    [Xt,Yt]=meshgrid(Xrange,Yrange);
    Xt=Xt.*dxy;Yt=Yt.*dxy;
    Uq(i,:) = interp2(Xt,Yt,vel_u(min(Xrange)+1:max(Xrange)+1,min(Yrange)+1:max(Yrange)+1),Xq,Yq);
    Vq(i,:) = interp2(Xt,Yt,vel_v(min(Xrange)+1:max(Xrange)+1,min(Yrange)+1:max(Yrange)+1),Xq,Yq);
end

Initial_con=[ID,T,X,Y,Uq,Vq];

% moving and changing 

tmax=length(t);   

for i=1:n
    DT(i,1)=Initial_con(i,2);
    DX(i,1)=Initial_con(i,3);
    DY(i,1)=Initial_con(i,4);
    DU(i,1)=Initial_con(i,5);
    DV(i,1)=Initial_con(i,6);
    for ti=2:tmax
          
        DX(i,ti)=DX(i,ti-1)+DU(i,ti-1)*ts;
        DY(i,ti)=DY(i,ti-1)+DV(i,ti-1)*ts;
        DT(i,ti)=DT(i,ti-1)+ts;
        
        vel_u1=U{ti-1}'; 
        vel_v1=V{ti-1}'; 
        vel1=sqrt(vel_u1.^2+vel_v1.^2);

        vel_u2=U{ti}'; 
        vel_v2=V{ti}'; 
        vel2=sqrt(vel_u2.^2+vel_v2.^2);

        % vel interp
%         [DU(i,t), DV(i,t)]=SpatialTimeInterp(DT(i,t),DX(i,t),DY(i,t),vel_u1,vel_v1,vel_u2,vel_v2,tt_num(tIdx-1:tIdx));
%         [Uq Vq] = SpatialTimeInterp (DT,DX,DY,vel_u1,vel_v1,vel_u2,vel_v2,tt_num)
%         [DU(i,ti), DV(i,ti)]=SpatialTimeInterp(DT(i,ti),DX(i,ti),DY(i,ti),vel_u,vel_v,vel_u,vel_v,tt_num(tIdx-1:tIdx));
%         nearby 4 points
        Xrange=[floor(DX(i,ti)/dxy):ceil(DX(i,ti)/dxy)];
        Yrange=[floor(DY(i,ti)/dxy):ceil(DY(i,ti)/dxy)];
        if max(Xrange)<length(x) & max(Yrange)<length(y) & min(Xrange)>=0 & min(Yrange)>=0
        Xq=DX(i,ti);Yq=DY(i,ti);
        [Xt,Yt]=meshgrid(Xrange,Yrange);
        Xt=Xt.*dxy;Yt=Yt.*dxy;
        DU(i,ti) = interp2(Xt,Yt,vel_u2(min(Xrange)+1:max(Xrange)+1,min(Yrange)+1:max(Yrange)+1),Xq,Yq);
        DV(i,ti) = interp2(Xt,Yt,vel_v2(min(Xrange)+1:max(Xrange)+1,min(Yrange)+1:max(Yrange)+1),Xq,Yq);
        else
        DU(i,ti)=0;
        DV(i,ti)=0;
        end
    end

%     DD_lon(i,t)
%     DD_lat(i,t)
end

% plot the trajectories
figure(2)
for i=1:n
    plot(DX(i,:),DY(i,:)); hold on;
end
axis([0,6.29,0,3.14])
title('Particle trajectories')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Below is the GPR predicted velocity fields %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch Casenumber
    case 1
        % Below is dataset1 for a longer time interval 0.5
        load('DoubleGyre_50drifter_40timestep_DD_OI.mat');
        ts=0.5;
    case 2
        % Below is dataset2 for a shorter time interval 0.2
        load('DoubleGyre_50drifter_100timestep_DD_OI.mat');
        ts=0.2;   
end

yr = double(yr);
mon = double(mon);
hr = double(hr);
day = double(day);
mn = double(mn);

U0 = U;
V0 = V;
Ku0 = Ku;
Kv0 = Kv;
xg0 = xg0;
yg0 = yg0;

yob0 = yob;
xob0 = xob;

lat = latt;
lon = lont;
Vel=sqrt(U0.^2+V0.^2);

DateVector0 = [yr; mon; day; hr; mn; zeros(1,length(mn)) ]';
time0 = datestr(DateVector0);
[l_time0,~]=size(time0);
dxy=abs(xg0(1,1,2)-xg0(1,1,1));
skipX = 0.2/dxy; %spacing for arrows
skipY = 0.2/dxy; %spacing for arrows



for i =1 : length(tg)
    figure(3)
    clf(3)
%     title(['GPR predicted TimeStep:',num2str(tg(i))],'FontSize',10)
    axis([0,6.29,0,3.14]);
    box on
    min_xg0=min(min(xg0(i,:)));
    max_xg0=max(max(xg0(i,:)));
    min_yg0=min(min(yg0(i,:)));
    max_yg0=max(max(yg0(i,:)));
    set(gca,'FontSize',10)
    [I J]=size(squeeze(U(i,:,:)));
	hold on
    l=streamslice(squeeze(xg0(i,1:skipX:end,1:skipY:end)),squeeze(yg0(i,1:skipX:end,1:skipY:end)),squeeze(U0(i,1:skipX:end,1:skipY:end)),squeeze(V0(i,1:skipX:end,1:skipY:end)));
    set(l,'color','b')
    plot(xob0(i,:),yob0(i,:),'ko','Linestyle','none','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
%     quiver(squeeze(xg0(i,1:skipX:end,1:skipY:end)),squeeze(yg0(i,1:skipX:end,1:skipY:end)),squeeze(U0(i,1:skipX:end,1:skipY:end)),squeeze(V0(i,1:skipX:end,1:skipY:end)),'color','k')
    
    hold off
    pause(.3)
end
