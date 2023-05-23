clear; clc;
addpath(genpath('C:\Users\samxj\OneDrive - University of Miami\PhD submeso ML research\MATLAB\JFLAB'))
addpath(genpath('D:\PhD Research\ML-GPR'))

% GPR load
cd 'C:\Users\samxj\PycharmProjects\GPR_NCOM\Examples_from_Lodise'

GPRFile_name={'Test_1drifter_3h_DD_OI.mat'; 'Test_5drifter_3h_DD_OI.mat'; 'Test_10drifter_3h_DD_OI.mat';...
    'Test_20drifter_3h_DD_OI.mat'; 'Test_30drifter_3h_DD_OI.mat'; 'Test_40drifter_3h_DD_OI.mat';...
    'Test_50drifter_3h_DD_OI.mat'; 'Test_60drifter_3h_DD_OI.mat'; 'Test_70drifter_3h_DD_OI.mat';...
    'Test_80drifter_3h_DD_OI.mat'; 'Test_90drifter_3h_DD_OI.mat'; 'Test_100drifter_3h_DD_OI.mat';...
    'Vortex_150drifter_12h_DD_OI.mat';'Convergence_200drifter_6h_DD_OI.mat';
    };  
load(GPRFile_name{end})
% GPR
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

% plot
figure
hold on;
i = l_time0;
min_xg0=min(min(xg0(i,:)));
max_xg0=max(max(xg0(i,:)));
min_yg0=min(min(yg0(i,:)));
max_yg0=max(max(yg0(i,:)));
IndX=[floor(min_xg0)+1:ceil(max_xg0)+1];IndY=[floor(min_yg0)+1:ceil(max_yg0)+1];
set(gca,'FontSize',10)
s=pcolor(squeeze(xg0(i,1:skipX:end,1:skipY:end)),squeeze(yg0(i,1:skipX:end,1:skipY:end)),squeeze(Vel(i,1:skipX:end,1:skipY:end)));
axis equal;axis xy;

for i =1 : n
    plot(xob0(:,i),yob0(:,i),'-','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
    hold on;
end
i = l_time0;
plot(xob0(end,:),yob0(end,:),'or','MarkerSize',3,'MarkerEdgeColor','k','MarkerFacecolor','g')
% quiver(squeeze(xg0(i,1:skipX:end,1:skipY:end)),squeeze(yg0(i,1:skipX:end,1:skipY:end)),squeeze(U0(i,1:skipX:end,1:skipY:end)),squeeze(V0(i,1:skipX:end,1:skipY:end)),'color','k')
axis([floor(min_xg0) ceil(max_xg0) floor(min_yg0)-0.5 ceil(max_yg0)])

s.EdgeColor = 'none'; ; colorbar; caxis auto;
colormap parula
%     caxis([0 2]);

% title(['GPR ',time0(i,:)],'FontSize',10)
title(['GPR'],'FontSize',10)
box on
xlabel('X')
ylabel('Y')
hold off
axis_range=[ceil(min_xg0) floor(max_xg0) ceil(min_yg0) floor(max_yg0)];
axis equal;axis xy;
axis(axis_range)
caxis([0 2]);
%%
% NCOM plot
File_Dir=dir('D:\PhD Research Datasets\ML-GPR\Original dataset_NCOM_surface_velocity');
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
vel=sqrt(vel_u.^2+vel_v.^2); 
vel_ut=transpose(vel_u);vel_vt=transpose(vel_v);velt=transpose(vel);
X=[0:1799];Y=[0:1419];
[MX,MY]=meshgrid(X,Y);




