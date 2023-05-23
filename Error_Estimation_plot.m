% error estimation
% there is a problem for the velocity field plot, left bottom or center.
% fix later

clear; clc;close all;
cd 'C:\Users\samxj\OneDrive - University of Miami\PhD submeso ML research\MATLAB\GPR_NCOM\GPR_VelF_output'


% GPR data
load('Convergence_200drifter_6h_DD_OI_FIL.mat')
% NCOM data
File_Dir=dir('D:\PhD Research Datasets\ML-GPR\Original dataset_NCOM_surface_velocity\');
path=['D:\PhD Research Datasets\ML-GPR\Original dataset_NCOM_surface_velocity\'];

i=1;
for k=3:length(File_Dir)
    File_name(i)={File_Dir(k).name(1:end)};
    yy=str2num(File_Dir(k).name(10:13));
    mm=str2num(File_Dir(k).name(14:15));
    dd=str2num(File_Dir(k).name(16:17));
    hh=str2num(File_Dir(k).name(23:24));
    tt(i,:)=[yy mm dd hh];
    tt_num(i)=datenum(yy,mm,dd,hh,0,0);
    i=i+1;
end

FIdx=1;
NCOMlon=ncread([path,char(File_name(FIdx))],'longitude');
NCOMlat=ncread([path,char(File_name(FIdx))],'latitude'); 
vel_u=ncread([path,char(File_name(FIdx))],'vel_u'); 
vel_v=ncread([path,char(File_name(FIdx))],'vel_v'); 
vel_ut=transpose(vel_u);vel_vt=transpose(vel_v);
vel=sqrt(vel_u.^2+vel_v.^2); velt=transpose(vel);
X=[0:1799];Y=[0:1419];
[MX,MY]=meshgrid(X,Y);
[MLon,MLat]=meshgrid(NCOMlon-360,NCOMlat);

figure
s=pcolor(MLon,MLat,velt);
s.EdgeColor = 'none'; s.FaceColor = 'interp';
load coastlines
hold on
% plot(coastlon,coastlat,'r','linewidth',2)
geoshow(coastlat,coastlon,"DisplayType","polygon", ...
    "FaceColor",[.7 .7 .7])


s_x=[min(min(lont))+0.1,min(min(lont))+0.1,max(max(lont))+0.1,max(max(lont))+0.1,min(min(lont))+0.1];
s_y=[min(min(latt))+0.35,max(max(latt))+0.35,max(max(latt))+0.35,min(min(latt))+0.35,min(min(latt))+0.35];
plot(s_x,s_y,'r','linewidth',3)
% [min(lont),min(latt)]
% [min(lont),max(latt)]
% [max(lont),max(latt)]
% [max(lont),min(latt)]
% [min(lont),min(latt)]



% deal with the GPR data
[~,n]=size(xob);

yr = double(yr); mon = double(mon);
hr = double(hr); day = double(day);
mn = double(mn);
U0 = U; V0 = V;
Ku0 = Ku; Kv0 = Kv; 
xg0 = xg0; yg0 = yg0;
lat_degs = lat_degs; lon_degs = lon_degs;
yob0 = yob; xob0 = xob;
lat = latt; lon = lont;
Vel=sqrt(U0.^2+V0.^2);
DateVector0 = [yr; mon; day; hr; mn; zeros(1,length(mn)) ]';
time0 = datestr(DateVector0);
[l_time0,~]=size(time0);

skipX = 5; %spacing for arrows
skipY = 5; %spacing for arrows

i = l_time0;
min_xg0=min(min(xg0(i,:)));
max_xg0=max(max(xg0(i,:)));
min_yg0=min(min(yg0(i,:)));
max_yg0=max(max(yg0(i,:)));
IndX=[floor(min_xg0)+1:skipX:ceil(max_xg0)+1];IndY=[floor(min_yg0)+1:skipY:ceil(max_yg0)+1];

axis_range=[ceil(min_xg0) floor(max_xg0) ceil(min_yg0) floor(max_yg0)];
%%
% NCOM plot
figure
hold on
% imagesc(X,Y,transpose(vel));axis equal;axis xy;colorbar; 
s=pcolor(X(IndX),Y(IndY),velt(IndY,IndX));axis equal;axis xy;colorbar;
% quiver(MX(IndY,IndX),MY(IndY,IndX),vel_ut(IndY,IndX),vel_vt(IndY,IndX),'k'); 
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
plot(xob0(end,:),yob0(end,:),'or','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
axis equal;axis xy;
axis(axis_range)
title('NCOM','FontSize',10)
s.EdgeColor = 'none'; s.FaceColor = 'interp';

% GPR plot
figure
hold on;
set(gca,'FontSize',10)
i = l_time0;
s=pcolor(squeeze(xg0(i,1:skipX:end,1:skipY:end)),squeeze(yg0(i,1:skipX:end,1:skipY:end)),squeeze(Vel(i,1:skipX:end,1:skipY:end)));
axis equal;axis xy;

for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
i = l_time0;
plot(xob0(end,:),yob0(end,:),'or','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
% quiver(squeeze(xg0(i,1:skipX:end,1:skipY:end)),squeeze(yg0(i,1:skipX:end,1:skipY:end)),squeeze(U0(i,1:skipX:end,1:skipY:end)),squeeze(V0(i,1:skipX:end,1:skipY:end)),'color','k')
axis(axis_range)

s.EdgeColor = 'none';s.FaceColor = 'interp' ; colorbar; caxis auto;
colormap parula
%     caxis([0 2]);

% title(['GPR ',time0(i,:)],'FontSize',10)
title(['GPR'],'FontSize',10)
box on
xlabel('X')
ylabel('Y')
hold off
%     print(gcf,[time0(i,:) ],'-djpeg', '-r300' )
    
%% Ku plot
figure 
hold on
i = l_time0;
pcolor(squeeze(xg0(i,1:skipX:end,1:skipY:end)),squeeze(yg0(i,1:skipX:end,1:skipY:end)),squeeze(Ku(i,1:skipX:end,1:skipY:end)))
axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
axis(axis_range)
colorbar
colormap parula
% Kv plot
title(['Ku'],'FontSize',10)
box on
xlabel('X')
ylabel('Y')
hold off

%% Kv plot
figure 
hold on
i = l_time0;
pcolor(squeeze(xg0(i,1:skipX:end,1:skipY:end)),squeeze(yg0(i,1:skipX:end,1:skipY:end)),squeeze(Kv(i,1:skipX:end,1:skipY:end)))
axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
axis(axis_range)
colorbar
colormap parula
% Kv plot
title(['Kv'],'FontSize',10)
box on
xlabel('X')
ylabel('Y')
hold off

%% comparasion
tsc=6 % according to time0(tsc,:)
% FIdxc=2; % 6 hour according to File_name(FIdxc)
FIdxc=1; % 6 hour according to File_name(FIdxc)
vel_uc=ncread([path,char(File_name(FIdxc))],'vel_u'); 
vel_vc=ncread([path,char(File_name(FIdxc))],'vel_v'); 


% [GPR_X, GPR_Y]=latlon2local(squeeze(lat_degs(tsc,1:end,1:end)),squeeze(lon_degs(tsc,1:end,1:end)),0,origin);
% GPR_X=GPR_X./1000;
% GPR_Y=GPR_Y./1000;
GPR_X=squeeze(xg0(tsc,:,:));
GPR_Y=squeeze(yg0(tsc,:,:));
GPR_U=squeeze(U0(tsc,1:end,1:end));
GPR_V=squeeze(V0(tsc,1:end,1:end));
[I J]=size(GPR_U);
% becareful about rhe vel_uc here
for i=1:I
    for j=1:J
        Xrange=[floor(GPR_X(i,j))-1:ceil(GPR_X(i,j))+1];
        Yrange=[floor(GPR_Y(i,j))-1:ceil(GPR_Y(i,j))+1];
        [Xt,Yt]=meshgrid(Xrange,Yrange);
        NCOM_U(i,j)=interp2(Xt,Yt,vel_uc(min(Xrange)+1:max(Xrange)+1,min(Yrange)+1:max(Yrange)+1),GPR_X(i,j),GPR_Y(i,j));
        NCOM_V(i,j)=interp2(Xt,Yt,vel_vc(min(Xrange)+1:max(Xrange)+1,min(Yrange)+1:max(Yrange)+1),GPR_X(i,j),GPR_Y(i,j));
    end
end
%% absolute error
E_Diff_U=NCOM_U-GPR_U;
E_Diff_V=NCOM_V-GPR_V;
E_Diff_Vel=(sqrt(NCOM_U.^2+NCOM_V.^2)-sqrt(GPR_U.^2+GPR_V.^2));

% exact error u
figure
s=pcolor(GPR_X,GPR_Y,E_Diff_U);
set(s,'EdgeColor','none');
hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
% plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
title(["U Absolute Error (m/s) "])
colorbar; caxis auto

% exact error v
figure
s=pcolor(GPR_X,GPR_Y,E_Diff_V);
set(s,'EdgeColor','none');
hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
% plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
title(["V Absolute Error (m/s) "])
colorbar; caxis auto

% exact error vel
figure
s=pcolor(GPR_X,GPR_Y,E_Diff_Vel);
hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
% plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
title(["Vel Absolute Error (m/s) "])
colorbar; caxis auto

%% scatter plot Abs Err VS. ErrQ
Ku1=reshape(squeeze(Ku0(tsc,:,:)),[],[1]);
Kv1=reshape(squeeze(Kv0(tsc,:,:)),[],[1]);
figure
c=reshape(NCOM_U,[],[1]);
Idx = find(abs(c)>0.1);
Idx = Idx(1:2:length(Idx));
c=c(Idx);
r_E_Diff_U = reshape((E_Diff_U),[],[1]);
scatter(Ku1(Idx),abs(r_E_Diff_U(Idx)),5,c,'filled');
xlabel('ErrQ');
ylabel('|Absolute Error|');
colorbar;caxis([0,1])
figure
c=reshape(NCOM_V,[],[1]);
Idx = find(abs(c)>0.1);
Idx = Idx(1:2:length(Idx));
c=c(Idx);
r_E_Diff_V = reshape((E_Diff_V),[],[1]);
scatter(Kv1(Idx),abs(r_E_Diff_V(Idx)),5,c,'filled');
xlabel('ErrQ');
ylabel('|Absolute Error|');
colorbar;caxis([0,2])


%% 
c=reshape(NCOM_V,[],[1]);
Idx = find(abs(c)>1 & abs(c)<2);
c=c(Idx);
r_E_Diff_V = reshape((E_Diff_V),[],[1]);
scatter(Kv1(Idx),abs(r_E_Diff_V(Idx)),5,c,'filled');
xlabel('ErrQ');
ylabel('|Absolute Error|');
colorbar;caxis([0,2])


%% linear regression




%% relative error
Diff_U=(NCOM_U-GPR_U)./NCOM_U*100;
Diff_V=(NCOM_V-GPR_V)./NCOM_V*100;
Diff_Vel=(sqrt(NCOM_U.^2+NCOM_V.^2)-sqrt(GPR_U.^2+GPR_V.^2))./sqrt(NCOM_U.^2+NCOM_V.^2)*100;

% U percent error
figure
s=pcolor(GPR_X,GPR_Y,Diff_U);
set(s,'EdgeColor','none');
hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
% title(["U Percent Error (%) ",datestr(tt_num(FIdxc))])
title(["U Relative Error (%) "])
colorbar; caxis([-100,100])
% V percent error
figure
s=pcolor(GPR_X,GPR_Y,Diff_V);
set(s,'EdgeColor','none');
hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
% title(["V Percent Error (%) ",datestr(tt_num(FIdxc))])
title(["V Relative Error (%) "])
colorbar; caxis([-100,100])

figure
s=pcolor(GPR_X,GPR_Y,Diff_Vel);hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
% title(["Vel Percent Error (%) ",datestr(tt_num(FIdxc))])
title(["Vel Relative Error (%) "])
colorbar; caxis([-100,100])

%%
figure
scatter(Ku1,reshape((Diff_U),[],[1]));
ylim([-200,200])
xlabel('ErrQ')
ylabel('Relative Error')
figure
scatter(Kv1,reshape((Diff_V),[],[1]));
ylim([-200,200])
xlabel('ErrQ')
ylabel('Relative Error')
%% histogram
for n1 = 0.0085
figure
hold on;
t = tiledlayout(1,1);
ax1 = axes(t);
% n1=0.0088;
h=histogram(ax1,abs(E_Diff_V),'BinWidth',n1);
ylabel('Number of data')
xlabel('Absolute Error m/s')
ax2 = axes(t);
n2=0.0005;
h=histogram(ax2,Kv1,'BinWidth',n2);
h.FaceColor = [0 0.5 0.5];
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ylabel('Number of data')
xlabel('ErrQ')
xlim(ax1,[0,40*n1]);
xlim(ax2,[0,40*n2]);

ylim(ax1,[0,7000]);
ylim(ax2,[0,7000]);

lgd = legend([ax1.Children(1) ax2.Children(1)], 'V Absolute Error', 'V ErrQ');
end
%%

figure
hold on;
t = tiledlayout(1,1);
ax1 = axes(t);
n1=0.001;
h=histogram(ax1,abs(E_Diff_U),'BinWidth',n1);
ylabel('Number of data')
xlabel('Absolute Error m/s')
ax2 = axes(t);
n2=0.00001;
h=histogram(ax2,Ku1,'BinWidth',n2);
h.FaceColor = [0 0.5 0.5];
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ylabel('Number of data')
xlabel('ErrQ')
xlim(ax1,[0,40*n1]);
xlim(ax2,[0,40*n2]);

ylim(ax1,[0,8000]);
ylim(ax2,[0,8000]);

lgd = legend([ax1.Children(1) ax2.Children(1)], 'U Absolute Error', 'U ErrQ');



%%
figure
histogram(abs(E_Diff_V),'BinWidth',0.005);
hold on

histogram(Kv1,'BinWidth',0.0005)







%% magnitude error vs. exact/predeicted velocity
addpath 'C:\Users\samxj\OneDrive - University of Miami\PhD submeso ML research\MATLAB\GPR_NCOM\'
% U
Diff_U=(NCOM_U-GPR_U)./NCOM_U*100;
Diff_V=(NCOM_V-GPR_V)./NCOM_V*100;
Diff_Vel=(sqrt(NCOM_U.^2+NCOM_V.^2)-sqrt(GPR_U.^2+GPR_V.^2))./sqrt(NCOM_U.^2+NCOM_V.^2)*100;
% Exact
PlotPEvsVEL(NCOM_U, Diff_U)
title('U Relative Error vs. Absolute Velocity')
xlabel('Exact Velocity (log m/s)')

PlotPEvsVEL(NCOM_V, Diff_V)
title('V Relative Error vs. Absolute Velocity')
xlabel('Exact Velocity (log m/s)')
% Predicted
PlotPEvsVEL(GPR_U, Diff_U)
title('U Relative Error vs. Predicted Velocity')
xlabel('Predicted Velocity (log m/s)')

PlotPEvsVEL(GPR_V, Diff_V)
title('V Relative Error vs. Predicted Velocity')
xlabel('Predicted Velocity (log m/s)')
%%
% domain error (nearby 8 points) 
% i-1 to i+1 j-1 to j+1
E_Diff_U=NCOM_U-GPR_U;
E_Diff_V=NCOM_V-GPR_V;
P_Diff_U=(NCOM_U-GPR_U)./NCOM_U*100;
P_Diff_V=(NCOM_V-GPR_V)./NCOM_V*100;
Domain_E_Diff_U=E_Diff_U;
Domain_E_Diff_V=E_Diff_V;

[I J]=size(E_Diff_U);
for i=2:I-1
    for j=2:J-1
        Domain_E_Diff_U(i,j)=mean(mean(E_Diff_U(i-1:i+1,j-1:j+1)));
        Domain_E_Diff_V(i,j)=mean(mean(E_Diff_V(i-1:i+1,j-1:j+1)));
    end
end
        
Domain_P_Diff_U=Domain_E_Diff_U./NCOM_U*100;
Domain_P_Diff_V=Domain_E_Diff_V./NCOM_V*100;

% Exact error u
figure
s=pcolor(GPR_X,GPR_Y,Domain_E_Diff_U);
set(s,'EdgeColor','none');
hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
% plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
title(["Domain mean U Exact Error (m/s) "])
colorbar; caxis auto

% Exact error v
figure
s=pcolor(GPR_X,GPR_Y,Domain_E_Diff_V); 
set(s,'EdgeColor','none');
hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
% plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
title(["Domain mean V Exact Error (m/s) "])
colorbar; caxis auto


% U percent error
figure
s=pcolor(GPR_X,GPR_Y,Domain_P_Diff_U); hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
% title(["U Percent Error (%) ",datestr(tt_num(FIdxc))])
title(["Domain mean U Percent Error (%) "])
colorbar; caxis([0,20])
% V percent error
figure
s=pcolor(GPR_X,GPR_Y,Domain_P_Diff_V);hold on;axis equal;axis xy;
for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
plot(xob0(end,:),yob0(end,:),'or','MarkerSize',5)
axis(axis_range)
% title(["V Percent Error (%) ",datestr(tt_num(FIdxc))])
title(["Domain mean V Relative Error (%) "])
colorbar; caxis([-20,20])


%% density plot

% GPRFile_name={'Test_1drifter_3h_DD_OI.mat'; 'Test_5drifter_3h_DD_OI.mat'; 'Test_10drifter_3h_DD_OI.mat';...
%     'Test_20drifter_3h_DD_OI.mat'; 'Test_30drifter_3h_DD_OI.mat'; 'Test_40drifter_3h_DD_OI.mat';...
%     'Test_50drifter_3h_DD_OI.mat'; 'Test_60drifter_3h_DD_OI.mat'; 'Test_70drifter_3h_DD_OI.mat';...
%     'Test_80drifter_3h_DD_OI.mat'; 'Test_90drifter_3h_DD_OI.mat'; 'Test_100drifter_3h_DD_OI.mat'; };  

% GPRFile_name={'Convergence_20drifter_6h_DD_OI.mat';'Convergence_50drifter_6h_DD_OI.mat';'Convergence_100drifter_6h_DD_OI.mat';'Convergence_200drifter_6h_DD_OI.mat'};
GPRFile_name={'Convergence_200drifter_6h_DD_OI.mat'};
% GPRFile_name={'Vortex_20drifter_12h_DD_OI.mat';'Vortex_50drifter_12h_DD_OI.mat';'Vortex_100drifter_12h_DD_OI.mat';'Vortex_150drifter_12h_DD_OI.mat'};

for i_load=1:length(GPRFile_name)
    load(GPRFile_name{i_load})
        [~,n]=size(xob);

    yr = double(yr); mon = double(mon);
    hr = double(hr); day = double(day);
    mn = double(mn);
    U0 = U; V0 = V;
    Ku0 = Ku; Kv0 = Kv; 
    xg0 = xg0; yg0 = yg0;
    lat_degs = lat_degs; lon_degs = lon_degs;
    yob0 = yob; xob0 = xob;
    lat = latt; lon = lont;
    Vel=sqrt(U0.^2+V0.^2);
    DateVector0 = [yr; mon; day; hr; mn; zeros(1,length(mn)) ]';
    time0 = datestr(DateVector0);
    [l_time0,~]=size(time0);
    skipX = 1; %spacing for arrows
    skipY = 1; %spacing for arrows
    
    tsc=6 % according to time0(tsc,:)
    % FIdxc=2; % 6 hour according to File_name(FIdxc)
    FIdxc=1; % 6 hour according to File_name(FIdxc)
    vel_uc=ncread([path,char(File_name(FIdxc))],'vel_u'); 
    vel_vc=ncread([path,char(File_name(FIdxc))],'vel_v'); 


    % [GPR_X, GPR_Y]=latlon2local(squeeze(lat_degs(tsc,1:end,1:end)),squeeze(lon_degs(tsc,1:end,1:end)),0,origin);
    % GPR_X=GPR_X./1000;
    % GPR_Y=GPR_Y./1000;
    GPR_X=squeeze(xg0(tsc,:,:));
    GPR_Y=squeeze(yg0(tsc,:,:));
    GPR_U=squeeze(U0(tsc,1:end,1:end));
    GPR_V=squeeze(V0(tsc,1:end,1:end));
    [I J]=size(GPR_U);
    for i=1:I
        for j=1:J
            Xrange=[floor(GPR_X(i,j))-2:ceil(GPR_X(i,j))+2];
            Yrange=[floor(GPR_Y(i,j))-2:ceil(GPR_Y(i,j))+2];
            [Xt,Yt]=meshgrid(Xrange,Yrange);
            NCOM_U(i,j)=interp2(Xt,Yt,vel_uc(min(Xrange)+1:max(Xrange)+1,min(Yrange)+1:max(Yrange)+1),GPR_X(i,j),GPR_Y(i,j));
            NCOM_V(i,j)=interp2(Xt,Yt,vel_vc(min(Xrange)+1:max(Xrange)+1,min(Yrange)+1:max(Yrange)+1),GPR_X(i,j),GPR_Y(i,j));
        end
    end
    % percent error
    Diff_U=(NCOM_U-GPR_U)./NCOM_U*100;
    Diff_V=(NCOM_V-GPR_V)./NCOM_V*100;
    Diff_Vel=(sqrt(NCOM_U.^2+NCOM_V.^2)-sqrt(GPR_U.^2+GPR_V.^2))./sqrt(NCOM_U.^2+NCOM_V.^2)*100;

    [I J]=size(xob);
    n_sample_points(i_load)=I*J;
    Area(i_load)=(max(max(GPR_X))-min(min(GPR_X)))*(max(max(GPR_Y))-min(min(GPR_Y)));
    Density(i_load)=n_sample_points(i_load)/Area(i_load); % in /km^2
    Diff_Vel_inside(i_load)=mean(mean(abs(Diff_Vel(2:end-1,2:end-1))));
    clear NCOM_U NCOM_V
end
figure
plot(Density,Diff_Vel_inside,'o-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
xlabel('Density (n/km^2)');ylabel('Overall Relative Error (%)')