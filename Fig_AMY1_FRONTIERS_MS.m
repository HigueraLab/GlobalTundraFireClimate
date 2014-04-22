%% Fig_AMY1_FRONTIERS_MS
%% INITIALIZE WORKSPACE
clear all; close all; clc;
%%
gr = flipud(gray);
%%
% IMPORT SPATIAL DATA TO MAKE PANEL A (GLOBAL TUNDRA AND FIRE MAP)
cd L:\1_projectsData\Data_GIS\Arctic\CAVM\cavm_veg_shp\physiog_dd
cavm_veg = shaperead('physiog_dd.shp');   % cavm veg. shp file
cavm_veg_xy = [cavm_veg.X; cavm_veg.Y];   % cvam_veg polygons

% cd L:\3_labMembers\Adam\Proposals\2014\NESSF_2014\mapdata;
% terr = shaperead('terr_dd.shp'); % Terrestrial Landcover

% cd L:\3_labMembers\Adam\Proposals\2014\NESSF_2014\mapdata\official_teow\official
% biomes = shaperead('wwf_terr_ecos.shp'); % WWF Biome Classifications
% cavm_veg_xy = [cavm_veg.X; cavm_veg.Y];   % cavm_veg polygons

cd L:\1_projectsData\AK_tundra_project\AK_tundra_fire_JFSP_proposal\TundraFireJFSP_GIS\AK_tundra;
color = xlsread('cavm_lookup_table.xls','i2:k22');  % colors for CAVM 

cd L:\3_labMembers\Adam\Proposals\2014\NESSF_2014\mapdata\;
fire = double(geotiffread('fire_2001_2013_cavm_v2.tif'));
fireInfo = geotiffinfo('fire_2001_2013_cavm_v2.tif');

fire(fire==min(fire(:))) = NaN;
fire(fire==0) = NaN;

worldlo = importdata('L:\1_projectsData\Data_GIS\World_and_US\worldlo.mat');
% CREATE [LAT,LON] GRID FOR FIRE RASTER
[Lon,Lat] = meshgrid(-179.875:0.25:179.875,89.875:-0.25:-89.875);
% USER SET-UP PARAMETERS:
lat_lim = [40 90];      % latitude limit for map
long_lim = [-180 180];  % longitude limit for map
frame_on = {'on'};      % frame on = 1, off = 0;

lat_grid = [5];%[67.1];         % latitude locations for map grid
long_grid = [10];%[-153];       % longitude locations for map grid

site_lat = [66.8 67.4];         % latitude for site map
site_long = [-154.8 -150.8];    

% PANEL A: GLOBAL TUNDRA AND AREA BURNED
close all; figure (1); clf; set(gcf,'color','w','units','inches',...
    'position',[1 2 6.7 8])
subplot(4,4,[1 2 5 6]);
base_map = axesm('ortho','maplatlim',lat_lim,'maplonlim',long_lim,'frame',...
    char(frame_on));

% gridm ('GColor',[.3 .3 .3]); %tightmap
setm(base_map,'PLineLocation',[66.5],'PLabelMeridian',[120],...
    'Origin',[88.5 -121.5],'FFaceColor',[0.73 0.91 1.00]);%,'FLatLimit',([-Inf 40]))
   % coast line
% map = patchesm(lat,long,'edgecolor',[.8 .8 .8],'facecolor',[.8 .8 .8]); % coast line
% map = plotm(lat,long,'linewidth',1,'color',[.5 .5 .5]);   % coast line
hold on
geoshow(worldlo.POpatch,'facecolor','w','edgecolor',[.8 .8 .8])
zoom(2.1)
framem off
geoshow('landareas.shp', 'FaceColor', [.5 .5 .5],'edgecolor','none');
for i = 1:20
    cavm_veg_xy = [cavm_veg(i).X; cavm_veg(i).Y];
    col_in = cavm_veg(i).physiog;
    patchesm(cavm_veg_xy(2,:),cavm_veg_xy(1,:),'edgecolor','none',...
    'facecolor',color(col_in,:));
end
set(gca,'units','inches','position',[0.5 5 2.9 2.9]);
freezeColors;
hold on;
%
surfm(Lat,Lon,fire);
colormap([1 0 0]);
freezeColors;
hgrat = gridm('on'); 
set(hgrat,'Clipping','on')
textm(53.5,104,'a','FontName','Arial','FontWeight','Bold','FontSize',12);

%
% scaleruler ('units','km')   % Make map scale
% setm(handlem('scaleruler1'),'YLoc',-0.6,'XLoc',-0.50,'RulerStyle',...
%     'patches','MajorTick',[0 750 1500],'MinorTick',[0 250 500],...
%     'MajorTickLength',km2nm(300))
% PANEL B: ALASKA
% Load Data
a = 'L:\3_labMembers\Adam\Fire_Climate_Veg_Modeling\Data_Sources_Maps\2000_meter\TIF_FILES\VEG_Landscape_Features\AK_MASK_new.tif';
tifInfo = geotiffinfo(a);
[Lat,Lon] = geoGrd(tifInfo.RefMatrix,tifInfo.Height,tifInfo.Width,...
    tifInfo,0);

% Load Fire Data
cd L:\3_labMembers\Adam\Proposals\2014\NESSF_2014\mapdata\;
fire = shaperead('ak_firehx_dd.shp');
fire_poly = [fire.X; fire.Y];
% Load Climate Data
cd L:\3_labMembers\Adam\Fire_Climate_Veg_Modeling\Data_Sources_Maps\2000_meter\TIF_FILES\Climatologies\;
tjja = geotiffread('MeanJJATemp_1971_2000.tif');
pjja = geotiffread('TotJJAPcp_1971_2000.tif');

% Load Landscape and Climate Data
cd L:\3_labMembers\Adam\Fire_Climate_Veg_Modeling\Data_Sources_Maps\2000_meter\TIF_FILES\VEG_Landscape_Features\;
veg  = geotiffread('akVeg_newBF.tif');
veg(veg==5) = NaN;
latPlot = [58 90];
lonPlot = [-168 -141];

subplot(4,4,[3 4 7 8]);
axesm('tranmerc','MapLatLimit',latPlot,'MapLonLimit',lonPlot,...
   'Frame','off','Grid','off','MeridianLabel','off', ...
   'ParallelLabel','off','fontsize',10,'PLineLocation',[55:5:90],...
    'MLineLocation',[-140:-10:-168])

AK = shaperead('usastatehi', 'UseGeoCoords', true,...
            'Selector',{@(name) strcmpi(name,'Alaska'), 'Name'});
h = geoshow(AK, 'FaceColor', [0.5 0.5 0.5],'EdgeColor','none',...
    'LineWidth',1);
box on
% framem off
axis tight
sc = scaleruler('on');
setm(handlem('scaleruler1'),'YLoc',1.02,'XLoc',.062,...
    'MajorTick',0:100:200,...
    'MinorTick',0:50:100,...
    'MajorTickLabel',{'0','','200'},...
    'FontSize',8,...
    'FontName','Arial')
set(gca,'units','inches','position',[3.7 5 2.9 2.9]);
h1 = surfm(Lat,Lon,veg);  
% cavmCmap = [23 229.5 191.25; 
%             112 255 0; 
%             177 177 57;
%             255 255 125; 
%             53 94 59; 
%             255 76.5 76.5];
cavmCmap = [23 229.5 191.25; 
            112 255 0; 
            177 177 57;
            255 255 125];
cavmCmap = cavmCmap./255;
colormap(cavmCmap); freezeColors;
% geoshow(fire, 'FaceColor', 'r','edgecolor','k','LineWidth',0.5);
patchesm(fire_poly(2,:),fire_poly(1,:),'FaceColor',...
    [1 0 0],'EdgeColor','k','LineWidth',0.5);
textm(69.8,-171,'b','FontName','Arial','FontWeight','Bold','FontSize',12);
% hold on;
%% Climate Space Map For Global Tundra
%% CREATE CLIMATOLOGIES FOR GLOBAL SUMMER TEMPERATURE AND SUMMER
cd L:\1_projectsData\Data_climate\CRU\climatedata\;
tmp     = importdata('cru_ts3.21.2001.2010.tmp.dat');
pre     = importdata('cru_ts3.21.2001.2010.pre.dat');
tmp1112 = importdata('cru_ts3.21.2011.2012.tmp.dat');
pre1112 = importdata('cru_ts3.21.2011.2012.pre.dat');

tmp = [tmp; tmp1112];
pre = [pre; pre1112];

tmp_array = NaN(360,720,12);
pre_array = NaN(360,720,12);
l = 1;
for i = 1:12
    tmp_i = NaN(360,720,3);
    pre_i = NaN(360,720,3);
    k = 1;
    for j = 1:12
        if j == 6 || j == 7 || j == 8
            tmp_i(:,:,k) = tmp(l:(l+359),:);
            pre_i(:,:,k) = pre(l:(l+359),:);
            l = l + 360;
            k = k + 1;
        else 
            l = l + 360;
            continue
        end
    end
    tmp_array(:,:,i) = nanmean(tmp_i,3);
    pre_array(:,:,i) = nansum(pre_i,3);
end

tmp_avg = nanmean(tmp_array,3)./10;
tmp_avg(tmp_avg==-99.9) = NaN;

pre_avg = nanmean(pre_array,3)./10;
pre_avg(pre_avg==-299.700) = NaN;
%%
% subplot(2,2,3);

% Interpolate Climate Grids to same resolution as global fire grid
[Lon,Lat] = meshgrid(-179.75:0.5:179.75,89.75:-0.5:-89.75);
[LonI,LatI] = meshgrid(-179.875:0.25:179.875,89.875:-0.25:-89.875);
tmp_avg_I = interp2(Lon,Lat,tmp_avg,LonI,LatI,'nearest');
pre_avg_I = interp2(Lon,Lat,pre_avg,LonI,LatI,'nearest');

cd L:\3_labMembers\Adam\Proposals\2014\NESSF_2014\mapdata\;
veg = flipud(double(geotiffread('cavm_no_glacier.tif')));
veg(veg==max(veg(:))) = NaN;
veg(veg>0) = 1;
mask = veg;

cd L:\3_labMembers\Adam\Proposals\2014\NESSF_2014\mapdata\;
fire = flipud(double(geotiffread('fire_2001_2013_cavm_v2.tif')));
fire(fire==min(fire(:))) = NaN;

Z_X = tmp_avg_I.*mask; Z_Y = pre_avg_I.*mask;
fire = fire.*mask;
[Z_v_p, Z_v_f] = climateSpaceGrid(Z_X,Z_Y,mask,45,50,[0 12],[0 300],fire);


%% Initialize Climate Variables for X and Y axis
X  = Z_X;
Y  = Z_Y;
Z  = mask;

MinMax_x = [0 12];
MinMax_y = [0 300];
nXBin    = 45;
nYBin    = 50;
smooth   = 20;

%Vectorize variables
x = X(X>-9999); % Climate
y = Y(X>-9999); % Climate
z = Z(X>-9999); % Vegetation or Geographic Region

f = fire(X>-9999); % Fire or other Geographic Entity (occurrence,frequency,etc...)

[nYBin,nXBin] = size(Z_v_p);

minX = MinMax_x(1); maxX = MinMax_x(2);
bin_x = range(MinMax_x)/(nXBin-1);

minY = MinMax_y(1); maxY = MinMax_y(2);
bin_y = range(MinMax_y)/(nYBin-1);

xBin = minX:bin_x:maxX;
yBin = minY:bin_y:maxY;

% Smoothing Bin Size (interp2.m)
x_int = bin_x/smooth;
y_int = bin_y/smooth;
[XI YI] = meshgrid([minX:x_int:maxX],[minY:y_int:maxY]);

% Meshgrid for countour plot
[X_i Y_i]   = meshgrid(xBin,yBin);


% Linearly Interpolate Proportional Values of landscape
ZI_v = interp2(X_i,Y_i,Z_v_p,XI,YI,'linear');
%%
subplot(4,4,9);
f(f==0) = NaN;
boxplot([y,y.*f],'Symbol','w+',...
    'PlotStyle','traditional',...
    'Colors',[0.5 0.5 0.5; 1 0 0],...
    'Labels',{'' ''});
ylim([0 200]);
set(gca,'Units','Inches',...
    'Position',[0.5 2.6 0.6 2.3],...
    'YTickLabel',{'' '50' '100' '150' ''});
ylabel('Total Summer Precipitation (mm)',...
    'FontSize',10,'FontName','Arial');
%%
subplot(4,4,10);
%Contour fill - Landscape Proportion
[cF cHf] = contourf(XI,YI,ZI_v,10,'LineStyle','None');
% Set Colormap for filled contour plot
colormap(gr)
freezeColors;
xlim([1 13]);
ylim([0 200]);
set(gcf,'color',[1 1 1])
axis square;
% Smooth fire contours through 2-D interpolation
ZI_f = interp2(X_i,Y_i,Z_v_f,XI,YI,'linear');
hold on;
contBin = prctile(Z_v_f(Z_v_f>-9999),[1 100]);
contours = [contBin(1):range(contBin)/(6):contBin(2)];
cmap = [1.0 1.0 0.3;
        1.0 1.0 0.0;
        1.0 0.6 0.0;
        1.0 0.3 0.0;
        1.0 0.0 0.0;
        0.6 0.0 0.0;
        0.3 0.0 0.0];
for i = 1:length(contours)
    [mC hC] = contour(XI,YI,ZI_f,contours(i));
    set(hC,'Color',cmap(i,:),'LineWidth',1)
end
set(gca,'XTickLabel',{''},...
    'YTickLabel',{''},...
    'Units','inches',...
    'Position',[1.1 2.6 2.3 2.3]);
text(2,190,'c','FontName','Arial','FontSize',12,...
    'FontWeight','Bold');
%%
subplot(4,4,14);
boxplot([x,x.*f],'Symbol','w+','Orientation','horizontal',...
    'PlotStyle','traditional',...
    'Colors',[0.5 0.5 0.5; 1 0 0],...
    'Labels',{'' ''});
xlim([1 13]);
set(gca,'Units','Inches',...
    'Position',[1.1 2.0 2.3 0.6]);
xlabel('Mean Summer Temperature (\circC)',...
    'FontSize',10,'FontName','Arial');

%%
cd L:\3_labMembers\Adam\Fire_Climate_Veg_Modeling\FINAL_DATA_PROCESSING_AND_RESULTS\Data_Sources_Maps\Landscape_Veg\;
veg = (double(geotiffread('akVeg_newBF.tif')));
veg(veg==max(veg(:))) = NaN;
veg(veg>0) = 1;
mask = veg;

cd L:\3_labMembers\Adam\Fire_Climate_Veg_Modeling\FINAL_DATA_PROCESSING_AND_RESULTS\Data_Sources_Maps\Fire
fire = flipud(double(geotiffread('fire_1950_2009.tif')));

Z_X = tjja.*mask; Z_Y = pjja.*mask;
fire = fire.*mask;
[Z_v_p, Z_v_f] = climateSpaceGrid(Z_X,Z_Y,mask,45,50,[0 30],[0 700],fire);


% Initialize Climate Variables for X and Y axis
X  = Z_X;
Y  = Z_Y;
Z  = mask;

MinMax_x = [0 30];
MinMax_y = [0 700];
nXBin    = 45;
nYBin    = 50;
smooth   = 20;

%Vectorize variables
x = X(X>-9999); % Climate
y = Y(X>-9999); % Climate
z = Z(X>-9999); % Vegetation or Geographic Region

f = fire(X>-9999); % Fire or other Geographic Entity (occurrence,frequency,etc...)

[nYBin,nXBin] = size(Z_v_p);

minX = MinMax_x(1); maxX = MinMax_x(2);
bin_x = range(MinMax_x)/(nXBin-1);

minY = MinMax_y(1); maxY = MinMax_y(2);
bin_y = range(MinMax_y)/(nYBin-1);

xBin = minX:bin_x:maxX;
yBin = minY:bin_y:maxY;

% Smoothing Bin Size (interp2.m)
x_int = bin_x/smooth;
y_int = bin_y/smooth;
[XI YI] = meshgrid([minX:x_int:maxX],[minY:y_int:maxY]);

% Meshgrid for countour plot
[X_i Y_i]   = meshgrid(xBin,yBin);


% Linearly Interpolate Proportional Values of landscape
ZI_v = interp2(X_i,Y_i,Z_v_p,XI,YI,'linear');
%%
subplot(4,4,11);
f(f==0) = NaN;
boxplot([y,y.*f],'Symbol','w+',...
    'PlotStyle','traditional',...
    'Colors',[0.5 0.5 0.5; 1 0 0],...
    'Labels',{'' ''});
ylim([25 325]);
set(gca,'Units','Inches',...
    'Position',[3.7 2.6 0.6 2.3]);
%
subplot(4,4,12);
%Contour fill - Landscape Proportion
[cF cHf] = contourf(XI,YI,ZI_v,20,'LineStyle','None');
% Set Colormap for filled contour plot
colormap(gr)
freezeColors;
xlim([2 16]);
ylim([25 325]);
set(gcf,'color',[1 1 1])
axis square;
% Smooth fire contours through 2-D interpolation
ZI_f = interp2(X_i,Y_i,Z_v_f,XI,YI,'linear');
hold on;
contBin = prctile(Z_v_f(Z_v_f>-9999),[1 100]);
contours = [contBin(1):range(contBin)/(7):contBin(2)];
cmap = [1.0 1.0 0.6;
        1.0 1.0 0.3;
        1.0 1.0 0.0;
        1.0 0.6 0.0;
        1.0 0.3 0.0;
        1.0 0.0 0.0;
        0.6 0.0 0.0;
        0.3 0.0 0.0];
for i = 1:length(contours)
    [mC hC] = contour(XI,YI,ZI_f,contours(i));
    set(hC,'Color',cmap(i,:),'LineWidth',1)
end
set(gca,'XTickLabel',{''},...
    'YTickLabel',{''},...
    'Units','inches',...
    'Position',[4.3 2.6 2.3 2.3]);
text(3,310,'d','FontName','Arial','FontSize',12,...
    'FontWeight','Bold');

subplot(4,4,16);
boxplot([x,x.*f],'Symbol','w+','Orientation','horizontal',...
    'PlotStyle','traditional',...
    'Colors',[0.5 0.5 0.5; 1 0 0],...
    'Labels',{'' ''});
xlim([2 16]);

set(gca,'Units','Inches',...
    'Position',[4.3 2.0 2.3 0.6]);
xlabel('Mean Summer Temperature (\circC)',...
    'FontSize',10,'FontName','Arial');
