%This script will read in Supersucker gridded data for 2017 surveys in the
%Chukchi Sea, described in Juranek et al. (2023), doi: 10.1029/2022JC019292
%this script is designed to work with data downloaded to a local drive
%--------------------------------------------------------------------------
%TO USE:: update desired transect for plotting and scaling factor in 
% lines 121-124
%--------------------------------------------------------------------------
%CALLS: perceptually uniform colormaps in function cmocean, available from 
%Github (https://github.com/chadagreene/cmocean)
%Thyng, Kristen, et al. “True Colors of Oceanography: Guidelines for 
%Effective and Accurate Colormap Selection.” Oceanography, vol. 29, no. 3, 
%The Oceanography Society, Sept. 2016, pp. 9–13, doi:10.5670/oceanog.2016.66.
%--------------------------------------------------------------------------
%last updated Jan 28, 2023 by L. Juranek (laurie.juranek@oregonstate.edu)
%--------------------------------------------------------------------------
clear all
close all


%read in Supersucker files in netcdf format
fid='BC1proc.nc';
finfo=ncinfo(fid);
% loop and read in all variables
for j=1:numel(finfo.Variables)
    % extract the jth variable (type = string)
    var = finfo.Variables(j).Name;

    % use dynamic field name to add this to the structure
    BC1.(var) = ncread(fid,var);

    % convert from single to double, if that matters to you (it does to me)
    if isa(BC1.(var),'single')
        BC1.(var) = double(BC1.(var));
    end
end

fid='BC2proc.nc';
finfo=ncinfo(fid);
% loop over the variables
for j=1:numel(finfo.Variables)
    % extract the jth variable (type = string)
    var = finfo.Variables(j).Name;

    % use dynamic field name to add this to the structure
    BC2.(var) = ncread(fid,var);

    % convert from single to double, if that matters to you (it does to me)
    if isa(BC2.(var),'single')
        BC2.(var) = double(BC2.(var));
    end
end

fid='WT1proc.nc';
finfo=ncinfo(fid);
% loop over the variables
for j=1:numel(finfo.Variables)
    % extract the jth variable (type = string)
    var = finfo.Variables(j).Name;

    % use dynamic field name to add this to the structure
    WT1.(var) = ncread(fid,var);

    % convert from single to double, if that matters to you (it does to me)
    if isa(WT1.(var),'single')
        WT1.(var) = double(WT1.(var));
    end
end

fid='WT2proc.nc';
finfo=ncinfo(fid);
% loop over the variables
for j=1:numel(finfo.Variables)
    % extract the jth variable (type = string)
    var = finfo.Variables(j).Name;

    % use dynamic field name to add this to the structure
    WT2.(var) = ncread(fid,var);

    % convert from single to double, if that matters to you (it does to me)
    if isa(WT2.(var),'single')
        WT2.(var) = double(WT2.(var));
    end
end

fid='BHSproc.nc';
finfo=ncinfo(fid);
% loop over the variables
for j=1:numel(finfo.Variables)
    % extract the jth variable (type = string)
    var = finfo.Variables(j).Name;

    % use dynamic field name to add this to the structure
    BHS.(var) = ncread(fid,var);

    % convert from single to double, if that matters to you (it does to me)
    if isa(BHS.(var),'single')
        BHS.(var) = double(BHS.(var));
    end
end

fid='DBO4proc.nc';
finfo=ncinfo(fid);
% loop over the variables
for j=1:numel(finfo.Variables)
    % extract the jth variable (type = string)
    var = finfo.Variables(j).Name;

    % use dynamic field name to add this to the structure
    DBO.(var) = ncread(fid,var);

    % convert from single to double, if that matters to you (it does to me)
    if isa(DBO.(var),'single')
        DBO.(var) = double(DBO.(var));
    end
end

%this script is only designed to plot select variables for a single
%transect. To change transect uncomment relevant line below and comment
%others
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trans=BC1;printid='BC1';fac=.4; %choose tranect to plot, fac used to scale aspect of shorter BC transect
% trans=BHS;printid='BHS';fac=1;
% trans=DBO;printid='DBO';fac=1;
% trans=WT2;printid='WT2';fac=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate some things for plotting first
trans.dist=trans.distance; %this variable name is inconsistent among transects

%meshgrid depth and distance for pcolor plots
[dpthgrid,distgrid]=meshgrid(trans.depth, trans.dist);

%find N2 max -- lots of noise between 5-10 m so this eliminates this depth interval
[stab_max,sbmax_i]=max(trans.N2(:,10:end), [],2);
pycno=filloutliers(trans.depth(sbmax_i+10),'linear'); %the +10 is to align the N2 max index with the full depth grid (since 1:10 excluded above)
        
%set up parameters for bottom patch
yl=120; %ylimit for plot
bottom=trans.bottom_depth;
bottompoly=[bottom' yl yl bottom(1)];
bottompolx=[trans.dist' trans.dist(end)  trans.dist(1) trans.dist(1)];

f1=figure; %O2 sat
pcolor(distgrid,dpthgrid, 100*(trans.O2_sat-1))
shading flat
hold on
[c,h]=contour(distgrid,dpthgrid, trans.pdens-1000,[23.3:.4:27.7], 'color', .7*[1 1 1]);
plot(trans.dist,pycno,'k.')

ylabel('Depth')
xlabel('Transect Distance (km)')
set(gca, 'fontsize', 11, 'fontweight', 'bold')
set(gcf,'PaperPosition',[0.75 1 6 3])
cmocean('delta')
colorbar
caxis([-30 30])
ylim([0 yl])
hh=patch(bottompolx,bottompoly,'k');
set(hh,'facecolor',[.9 .9 .9],'edgecolor',[0 0 0])
set(gca,'layer','top')    
set(gca, 'fontsize', 11, 'fontweight', 'bold', 'ydir', 'reverse')
    set(gcf,'PaperPosition',[0.75 1 fac*6 3])
    
f2=figure; %N+N note NO3 refers to NO3 + NO2 > the analysis used cannot distinguish between these
pcolor(distgrid,dpthgrid, trans.NO3)
shading flat
hold on
[c,h]=contour(distgrid,dpthgrid, trans.pdens-1000,[23.3:.4:27.7], 'color', .7*[1 1 1]);
plot(trans.dist,pycno,'k.')
ylabel('Depth')
xlabel('Transect Distance (km)')
cmocean('matter')
colorbar
caxis([0 10])
ylim([0 yl])
hh=patch(bottompolx,bottompoly,'k');
set(hh,'facecolor',[.9 .9 .9],'edgecolor',[0 0 0])
set(gca,'layer','top')
set(gca, 'fontsize', 11, 'fontweight', 'bold', 'ydir', 'reverse')
    set(gcf,'PaperPosition',[0.75 1 fac*6 3])
    
f3=figure; % NH4
pcolor(distgrid,dpthgrid, trans.NH4)
shading flat
hold on
[c,h]=contour(distgrid,dpthgrid, trans.pdens-1000,[23.3:.4:27.7], 'color', .7*[1 1 1]);
plot(trans.dist,pycno,'k.')
ylabel('Depth')
xlabel('Transect Distance (km)')
cmocean('matter')
colorbar
caxis([0 6])
ylim([0 yl])
hh=patch(bottompolx,bottompoly,'k');
set(hh,'facecolor',[.9 .9 .9],'edgecolor',[0 0 0])
set(gca,'layer','top')
set(gca, 'fontsize', 11, 'fontweight', 'bold', 'ydir', 'reverse')
    set(gcf,'PaperPosition',[0.75 1 fac*6 3])
    

f4=figure; %beam transmission
pcolor(distgrid,dpthgrid, trans.beamc)
shading flat
hold on
[c,h]=contour(distgrid,dpthgrid, trans.pdens-1000,[23.3:.4:27.7], 'color', .7*[1 1 1]);
plot(trans.dist,pycno,'k.')
ylabel('Depth')
xlabel('Transect Distance (km)')
cmocean('-rain')
colorbar
caxis([0 1])
ylim([0 yl])
hh=patch(bottompolx,bottompoly,'k');
set(hh,'facecolor',[.9 .9 .9],'edgecolor',[0 0 0])
set(gca,'layer','top')
set(gca, 'fontsize', 11, 'fontweight', 'bold', 'ydir', 'reverse')
    set(gcf,'PaperPosition',[0.75 1 fac*6 3])
    


f5=figure; %chl voltage
%note that the fluorescence sensor was not rigorously calibrated and
%experienced drift issues throughout the cruise. We recommend that these data only be used in a
%relative or qualitative sense
%set dark value based on data qc:
if printid=='BHS' %dark value drifts during BHS transect
    dark=min(trans.Chl_voltage,[],2);
elseif printid=='WT2'
    dark=0.10;
elseif printid=='DBO' %this transect had issues with intermittent spikes in chl voltage that were not real
    dark=0;
    j=nansum(trans.Chl_voltage,2);
    id=find(j>3.5); %removing spikes in data
    trans.Chl_voltage(id,:)=NaN;
else
    dark=0;
end

pcolor(distgrid,dpthgrid, log10(abs(trans.Chl_voltage-dark)))
shading flat
hold on
[c,h]=contour(distgrid,dpthgrid, trans.pdens-1000,[23.3:.4:27.7], 'color', .7*[1 1 1]);
plot(trans.dist,pycno,'k.')
ylabel('Depth')
xlabel('Transect Distance (km)')
cmocean('-deep')
c=colorbar;
caxis(log10([0.05 0.5]))
ylim([0 yl])
hh=patch(bottompolx,bottompoly,'k');
set(hh,'facecolor',[.9 .9 .9],'edgecolor',[0 0 0])
set(gca,'layer','top')
set(gca, 'fontsize', 11, 'fontweight', 'bold', 'ydir', 'reverse')
    set(gcf,'PaperPosition',[0.75 1 fac*6 3])
    
f6=figure; %salinity
pcolor(distgrid,dpthgrid, trans.salinity)
shading flat
hold on
[c,h]=contour(distgrid,dpthgrid, trans.pdens-1000,[23.3:.4:27.7], 'color', .7*[1 1 1]);
clabel(c,h);
plot(trans.dist,pycno,'k.')
ylabel('Depth')
xlabel('Transect Distance (km)')
cmocean('haline')
colorbar
caxis([29 33])
ylim([0 yl])
hh=patch(bottompolx,bottompoly,'k');
set(hh,'facecolor',[.9 .9 .9],'edgecolor',[0 0 0])
set(gca,'layer','top')
set(gca, 'fontsize', 11, 'fontweight', 'bold', 'ydir', 'reverse')
    set(gcf,'PaperPosition',[0.75 1 fac*6 3])
    
f7=figure; %temperature
pcolor(distgrid,dpthgrid, trans.temperature)
shading flat
hold on
[c,h]=contour(distgrid,dpthgrid, trans.pdens-1000,[23.3:.4:27.7], 'color', .7*[1 1 1]);
clabel(c,h);
plot(trans.dist,pycno,'k.')
ylabel('Depth')
xlabel('Transect Distance (km)')
cmocean('thermal')
colorbar
caxis([-2 10])
ylim([0 yl])
hh=patch(bottompolx,bottompoly,'k');
set(hh,'facecolor',[.9 .9 .9],'edgecolor',[0 0 0])
set(gca,'layer','top')
set(gca, 'fontsize', 11, 'fontweight', 'bold', 'ydir', 'reverse')
    set(gcf,'PaperPosition',[0.75 1 fac*6 3])
    
    
f8=figure; %turbulence
pcolor(distgrid,dpthgrid, log10(trans.KT))
shading flat
hold on
[c,h]=contour(distgrid,dpthgrid, trans.pdens-1000,[23.3:.4:27.7], 'color', .7*[1 1 1]);
clabel(c,h);
plot(trans.dist,pycno,'k.')
ylabel('Depth')
xlabel('Transect Distance (km)')
cmocean('thernal')
colorbar
caxis([-7 -2])
ylim([0 yl])
hh=patch(bottompolx,bottompoly,'k');
set(hh,'facecolor',[.9 .9 .9],'edgecolor',[0 0 0])
set(gca,'layer','top')
set(gca, 'fontsize', 11, 'fontweight', 'bold', 'ydir', 'reverse')
    set(gcf,'PaperPosition',[0.75 1 fac*6 3])
    

print('-f1',[printid,'_O2sat'],'-depsc');
print('-f2',[printid,'_nitrate'],'-depsc');
print('-f3',[printid,'_nh4'],'-depsc');
print('-f4',[printid,'_beamc'],'-depsc');
print('-f5',[printid,'_chlvolt'],'-depsc');
print('-f6',[printid,'_salinity_pdens'],'-depsc');
print('-f7',[printid,'_temp_pdens'],'-depsc');
print('-f8',[printid,'_KT'],'-depsc');

