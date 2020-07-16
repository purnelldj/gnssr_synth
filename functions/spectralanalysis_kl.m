function [t_rh,rh_adj,rh_nadj,rms_adj,rms_nadj] = spectralanalysis_kl(startdate,enddate,station,slvlrdir,tgstring,...
    redconstits,doelvlims,removeoutliers,makefig)

% this function takes spectral analysis estimates of SNR data and other
% stats from observations (obtained using analyzesnr_fun.m)
% and creates a time series of adjusted reflector heights following the
% method outlined in Larson et al. (2017) "A 10 year comparison..."

% written by David Purnell (2020)

% INPUTS
% startdate: in datenum format
% enddate: in datenum format
% station: string e.g., 'sc02'
% slvlrdir: path to directory containing observed or synthetic slvlr
% arrays (i.e., using analyzesnr_fun.m or makesnr_fun.m)
% tgstring: path to tide gauge data in format 'xaxis' 'slvl'
% if no tg, leave empty: ''
% redoncstits: if set to 1 then use full 145 tidal constituents for
% adjustment, otherwise just use 5, can choose which 5 below
% only set to 1 if time series is long (1 year +)
% doelvlims: set 1 to only use data within set limits of higher and lower elevation
% angle defined in station input file (recommended)
% removeoutliers: set 1 to remove outliers greater than 3 standard deviations away
% from smoothed time series (recommended)
% makefig: set 1 to make a figure

% OUTPUTS
% t_rh: time index of output reflector heights
% rh_adj: adjusted reflector heights
% rh_nadj: unadjusted reflector heights
% rms_adj: RMS of tide gauge and adjusted reflector heights (if tide gauge given)
% rms_nadj: RMS of tide gauge and unadjusted reflector heights (if tide gauge given)

% other parameters
l1=1; % L1 signal choose (1) or (0)
l2=0; % L2 signal choose (1) or (0)
% constellations
glo=0; % choose (1) or (0)
gps=1; % choose (1) or (0)
% for looking at the direction of satellite overpassess
justfwd=0; % 1 for fwd, 2 for not fwd, 0 for off (all data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pwdstr=pwd;
addpath([pwdstr,'/functions/'])
addpath([pwdstr,'/functions/t_tide_v1.4beta'])
run([pwdstr,'/functions/station_inputs/',station,'_input'])

% HERE IS WHERE YOU ARE LOADING ALL DATA TO BE ANALYZED
slvlrt=[];
datevect=startdate-1;
while datevect<enddate
    datevect=datevect+1;
if exist([slvlrdir,'/',num2str(datevect),'.mat'],'file')==2
    load([slvlrdir,'/',num2str(datevect),'.mat'])
    slvlrt=[slvlrt;slvlr];
    clear slvlr
end
end
slvlr=slvlrt;
clear slvlrt

if size(slvlr,1)<1
    disp('apparently there is no data in this time range')
    return
end

% to just look at satellites with increasing or decreasing elevation w/ t
if justfwd==1
out=slvlr(:,3)<0;
slvlr(out,:)=[];
elseif justfwd==2
out=slvlr(:,3)>0;
slvlr(out,:)=[];
end

% getting desired satellite constellations
% need to add galileo and beidou
if glo==0
    delete=slvlr(:,2)>32;
    slvlr(delete,:)=[];
end
if gps==0
    delete=slvlr(:,2)<33;
    slvlr(delete,:)=[];
end

% getting rid of bad elv lims
if doelvlims==1
    delete=abs(slvlr(:,4)-elv_low)>elvlims;
    slvlr(delete,:)=[];
    delete=abs(slvlr(:,5)-elv_high)>elvlims;
    slvlr(delete,:)=[];
end

% removing mean from time series
tempmean=nanmean(slvlr(:,7));
slvlr(:,7)=tempmean-slvlr(:,7);
tempmean=nanmean(slvlr(:,11));
slvlr(:,11)=tempmean-slvlr(:,11);

% GET RID OF NANS
if l1==1 && l2==1
keep=isnan(slvlr(:,11))==0;
tsize=size(slvlr,1)+1;
slvlr=[slvlr;slvlr(keep,:)];
slvlr(tsize:end,7)=slvlr(tsize:end,11);
slvlr(:,11)=NaN;
elseif l1==0 && l2==1
keep=isnan(slvlr(:,11))==0;
slvlr=slvlr(keep,:);
slvlr(:,7)=slvlr(:,11);
end
delete=isnan(slvlr(:,7))==1;
slvlr(delete,:)=[];
slvlr=sortrows(slvlr,1);


% NOW START DOING LEAST SQUARES TIDAL FIT
rh_nadj=slvlr(:,7); % sea level measurements
tanthter=slvlr(:,3)./3600; % to convert to per hour
t=slvlr(:,1);
% getting rid of outliers
if removeoutliers==1
hsmooth=smoothdata(rh_nadj,'movmean',5);
diff1=abs(rh_nadj-hsmooth);
std1=std(diff1);
delete=diff1(:,1)>3*std1; % here choose sigma bounds to remove
rh_nadj(delete,:)=[];
t(delete,:)=[];
tanthter(delete,:)=[];
end
pointsperday=numel(t)/(enddate+1-startdate);
disp(['points per day =',num2str(pointsperday)])
load('tidefreqs.mat')
if redconstits==0
ju=1:145;
else
ju=[12 20 41 47 56]; % O1, K1, N2, M2, S2
end
coefs_0=rand(numel(ju)*2,1)*2*sqrt(0.005)-sqrt(0.005);
freqs=freqs(ju);
names=names(ju,:);
tempfun=@(coefs) tidemod_kl(coefs,t,rh_nadj,ju,tanthter,freqs,plat);
options=optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective',...
    'Display','off'); % off for now
coefs_ls=lsqnonlin(tempfun,coefs_0,[],[],options); % here is the least squares
rh_adj=tidemod_kl_plot(coefs_ls,t,rh_nadj,ju,tanthter,freqs,plat);
rh_adj=rh_adj-nanmean(rh_adj);
rh_nadj=rh_nadj-nanmean(rh_nadj);
disp('summary of tides')
[A,G,names] = coef2ampphase(coefs_ls,names);
t_rh=t;

if numel(tgstring)>0
load(tgstring)
tidexp=xaxis;
tideyp=slvl;
tiderms=interp1(xaxis,slvl,t_rh,'linear');
tempmean=nanmean(tiderms);
tiderms=tiderms-tempmean;
tideyp=tideyp-tempmean;

tmps=isnan(rh_adj)==0 & isnan(tiderms)==0;
rms_adj=rms(tiderms(tmps)-rh_adj(tmps));
rms_nadj=rms(tiderms(tmps)-rh_nadj(tmps));

disp(['rms_adj=',num2str(rms_adj*100),' cm'])
disp(['rms_nadj=',num2str(rms_nadj*100),' cm'])
else
    rms_adj=NaN;
    rms_nadj=NaN;
end

if makefig==1
    %%
% PLOTTING
% Defaults
width =10;     % Width in inches % 3.5 was for putting two on one line i think
height =5;    % Height in inches
lw = 1.5;      % LineWidth
msz = 5;       % MarkerSize
% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*70]);
% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2; %
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
close all  

figure('visible','on')
plot(tidexp,tideyp,'k','linewidth',1.5)
hold on
scatter(t_rh,rh_nadj,'r','linewidth',1)
hold on
scatter(t_rh,rh_adj,'b','linewidth',1.5)
h=legend('tide gauge','LSP estimates','with correction');
axis([startdate enddate+1 -2 2.5])
fs=15;
datetick('x',1,'keeplimits','keepticks')
ylabel('Sea level (m)','interpreter','latex','fontsize',fs)
set(gca,'ticklabelinterpreter','latex','fontsize',fs)
set(h,'interpreter','latex','fontsize',fs)
print('spectralanalysis_testfig', '-dpng', '-r300');
end

end
