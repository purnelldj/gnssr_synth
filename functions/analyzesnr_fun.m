function [slvlr,lspy] = analyzesnr_fun(station,snr_data,tdatenum,satconsts,tropd,decimate,tempsnr,templsp,tempfresnel)

%%

% This function analyzes SNR data to get basic data for GNSS-Reflectometry
% sea level measurements
% Made by David Purnell (2020)
%
% INPUTS
%
% station: 4 letter character vector for identifier of GNSS station e.g., 'sc02'
% if you wish to add a new station, you need to add a new station input in
% the functions/station_inputs directory
% snr_data: array of snrdata in same format as output from rinex2snr.m
% tdatenum: day in matlab datenum format
% satconsts: 1 by 3 double set to 1 to include GPS, GLONASS, GALILEO
% E.g., [1 0 1] means GPS + GALILEO
% tropd: set to 1 for Williams et al. (2017) correction, 2 for
% Santamaria-Gomez et al. (2017)
% decimate: choose to decimate the data and analyse every [] seconds
% tempsnr: choose 1 to output figures of detrended synthetic SNR data to the
% directory 'tempoutput'
% templsp: choose 1 to output figures of Lomb-Scargle Periodogram data to the
% directory 'tempoutput'
% tempfresnel: choose 1 to output .kml figs of the start and end fresnel
% zones for each satellite arc to the directory 'tempoutput'
%
% OUTPUTS
% 
% slvlr: a (k,14) array of statistics from k satellite arcs
% column 1 is mean time of arc in datenum format
% column 2 is satellite # (GPS is 1-32, GLONASS 33-56, GALILEO 57-)
% column 3 is tan(elv)/d(elv)/dt term use for height corrections
% column 4 is min elv
% column 5 is max elv
% column 6 is mean azimuth
% column 7 is slvl measurement using L1 frequency band
% column 8 is the mean magnitude of SNR data in watt/watt
% column 9 is the peak of LSP
% column 10 is the variance of the SNR data
% column 11-14 is as in column 7-10 but for L2 frequency band
% 
% lspy: a (k,300) array of LSP spectrums for k satellite arcs
% the frequency is interpolated by default from 1:300, this can be changed
% by setting the value of lspy_flim
% 

% probably need to delete these
clear slvlr lspy

ltdrun=0; % change if you just want to produce a set number of figures

curdt=datetime(tdatenum,'convertfrom','datenum');
curjd=juliandate(curdt);
%strday=char(datetime(curdt,'format','DDD'));
%stryr=char(datetime(curdt,'format','yy'));

% choose satellite constellations to include
gps=satconsts(1);
glo=satconsts(2);
gal=satconsts(3);

% these should be your local directory (full strings)
pwdstr=pwd;
functionstr=[pwdstr,'/functions']; 
% directory that contains mpsim, lla2ecef, station_inputs 
% there is also a matlab function lla2ecef that is probably better

maxf1fix=0; % if you want to set a max freq and not worry about nyquist lim

if glo==1
    load([functionstr,'/glonasswlen.mat'])
end

addpath(functionstr)
addpath([functionstr,'/fresnel'])
run([functionstr,'/station_inputs/',station,'_input'])
if decimate~=0
    dt=decimate;
end
delete([pwdstr,'/tempoutput/*'])

[~,~,hell]=ecef2lla(staxyz(1),staxyz(2),staxyz(3));

cnt=0;
clear delete

% tropd
if tropd==1
curjd2=doy2jd(curyr,curday)-2400000.5;
hgtant=hell;
hgtlim(1)=hell-ahgt+ahgt_bounds;
hgtlim(2)=hell-ahgt-ahgt_bounds;
dlat(1) = plat*pi/180.d0;
dlon(1) = plon*pi/180.d0;
nstat = 1;
it = 0;
[pant,~,~,tmant,eant,ah,aw,lambda,~] = gpt2_1w (curjd2,dlat,dlon,hgtant,nstat,it);
[plim(1),~,~,tmlim(1),elim(1),~,~,~,~] = gpt2_1w (curjd2,dlat,dlon,hgtlim(1),nstat,it);
[plim(2),~,~,tmlim(2),elim(2),~,~,~,~] = gpt2_1w (curjd2,dlat,dlon,hgtlim(2),nstat,it);
end

%refraction
if tropd==2
curjd2=doy2jd(curyr,curday)-2400000.5;
hgtlim(1)=hell-ahgt-ahgt_bounds;
hgtlim(2)=hell-ahgt+ahgt_bounds;
dlat(1) = plat*pi/180.d0;
dlon(1) = plon*pi/180.d0;
nstat = 1;
it = 0;
[plim(1),tlim(1),~,~,~,~,~,~,~] = gpt2_1w (curjd2,dlat,dlon,hgtlim(1),nstat,it);
[plim(2),tlim(2),~,~,~,~,~,~,~] = gpt2_1w (curjd2,dlat,dlon,hgtlim(2),nstat,it);
end

% note SNR data is
% sat#, elev, azimuth, seconds, refl dist1, refl dist2, S1, S2
clear slvl

allsats=unique(snr_data(:,1)); 

for aa=1:size(allsats,1) %%%%%%%
    
    
    if gps==0 && allsats(aa)<33
        continue
    elseif glo==0 && (allsats(aa)>32 && allsats(aa)<57)
        continue
    elseif gal==0 && allsats(aa)>56
        continue
    end
    
ijk=0;
repgaps=1;
while ijk~=repgaps
    
ijk=ijk+1;
    
doo=allsats(aa);
tmpid=snr_data(:,1)==doo;
sat2=snr_data(tmpid,:);
in=sat2(:,3)>azi_low & sat2(:,3)<azi_high;
sat2trr=sat2(in,:);
in=sat2trr(:,2)<elv_high & sat2trr(:,2)>elv_low;
sat2tr=sat2trr(in,:);
sat2tr=sortrows(sat2tr,4);
if exist('azi_mask')==1
    out=sat2tr(:,3)>azi_mask(1) & sat2tr(:,3)<azi_mask(2);
    sat2tr(out)=[];
end
if decimate~=0
    modspl=mod(sat2tr(:,4),decimate);
    delt=modspl(:)>0;
    sat2tr(delt,:)=[];
end
%%%%%%%
times=sat2tr(:,4);
if size(times,1)<1
    continue
end
% to sort satellites which have more than one valid overpass
gaps=diff(times(:,1))>dt*10;
repgaps=sum(gaps)+1;
gapids=1:numel(gaps);
gaps=gapids(gaps);
if repgaps>1
    if ijk==1
        sat2tr(gaps(ijk)+1:end,:)=[];
    elseif ijk>1
        sat2tr(1:gaps(ijk-1),:)=[];
        if repgaps>2 && ijk<repgaps
            sat2tr(gaps(ijk)-gaps(ijk-1)+1:end,:)=[];
        end
    end
end
if size(sat2tr,1)<3
    continue
end

fwd=0;
if sat2tr(2,2)-sat2tr(1,2)>0
    fwd=1;
end
for a=3:size(sat2tr,1) % this is basically just deleting all data if it ...
    %gets to a point where it changes directions
    if sat2tr(a,2)-sat2tr(a-1,2)>0
        tmp=1;
    else
        tmp=0;
    end
    if tmp~=fwd
        times(a:end)=[];
        sat2tr(a:end,:)=[];
        break
    end
end
if fwd==0
    sat2tr=flipud(sat2tr);
end

if abs(sat2tr(end,4)-sat2tr(1,4))<300
    continue
end

nol1=0;
nol2=0;

sinelv=sind(sat2tr(:,2));

SNR1=sqrt(10.^(sat2tr(:,7)./10));
SNR1a=SNR1;
sinelv1=sinelv;
deln=isnan(SNR1(:,1))==1;
SNR1(deln,:)=[];
sinelv1(deln,:)=[];
p1=polyfit(sinelv1,SNR1,2);
y1=polyval(p1,sinelv);
SNR1dt=SNR1a-y1;
if allsats(aa)<33
L1car=(299792458/1575.42e06); % for GPS
elseif allsats(aa)<32+25
L1car=glonasswlen(allsats(aa)-32);
else
L1car=(299792458/1575.420e06); % for galileo
end
if maxf1fix==0
maxf1=numel(sinelv)/(2*(max(sinelv)-min(sinelv)));
else
maxf1=maxf1fix; 
end
prec1=0.001;
ovs=round(L1car/(2*prec1*(max(sinelv)-min(sinelv))));
fi=1:1:maxf1;
[psd,f,A2,~,~,~]=fLSPw(sinelv,SNR1dt,fi,0.05,ovs);
psd=cell2mat(psd);
f=cell2mat(f);
reflh1=f.*0.5*L1car;
[~,id]=max(psd(:));
pks=findpeaks(psd);
pks=sort(pks);

if tropd==1
% trop delay
preh1=reflh1(id);
hsfc=hell-preh1;
if hsfc<hgtlim(1) && hsfc>hgtlim(2)
psfc=interp1(hgtlim,plim,hsfc,'linear');
tmsfc=interp1(hgtlim,tmlim,hsfc,'linear');
esfc=interp1(hgtlim,elim,hsfc,'linear');
theta=sat2tr(:,2);
clear thetarefr
for jj=1:numel(theta)
tau=trop_delay_tp(curjd,plat,hell,hsfc,theta(jj),pant,tmant,eant,ah,aw,lambda,psfc,tmsfc,esfc);
thetarefr(jj)=asind(sind(theta(jj))+0.5*tau/preh1);
end
sinelv=sind(thetarefr).';
sinelv1=sinelv;
sinelv1(deln,:)=[];
p1=polyfit(sinelv1,SNR1,2);
y1=polyval(p1,sinelv);
SNR1dt=SNR1a-y1;
[psd,f]=plomb(SNR1dt,sinelv,maxf1,ovs,'normalized'); %
reflh1=f.*0.5*L1car;
[~,id]=max(psd(:));
pks=findpeaks(psd);
pks=sort(pks);
else
    thetarefr=sat2tr(:,2);
end

elseif tropd==2
% refraction
preh=reflh1(id);
hsfc=hell-preh;
psfc=interp1(hgtlim,plim,hsfc,'linear');
tsfc=interp1(hgtlim,tlim,hsfc,'linear');
theta=sat2tr(:,2);
dele = (1/60)*510*psfc/((9/5*tsfc+492)*1010.16).*cotd(theta+7.31./(theta+4.4));
sinelv=sind(theta+dele);
thetarefr=theta+dele;
sinelv1=sinelv;
sinelv1(deln,:)=[];
p1=polyfit(sinelv1,SNR1,2);
y1=polyval(p1,sinelv);
SNR1dt=SNR1a-y1;
[psd,f,A2,~,~,~]=fLSPw(sinelv,SNR1dt,fi,0.05,ovs);
psd=cell2mat(psd);
f=cell2mat(f);
reflh1=f.*0.5*L1car;
[~,id]=max(psd(:));

end
if reflh1(id) < ahgt-ahgt_bounds || reflh1(id) > ahgt+ahgt_bounds...
         || max(psd)<10*mean(pks(1:end-1))
    nol1=1;   
end

SNR2=sqrt(10.^(sat2tr(:,8)./10));
SNR2a=SNR2;
sinelv2=sinelv;
deln=isnan(SNR2(:,1))==1;
SNR2(deln,:)=[];
sinelv2(deln,:)=[];
if numel(sinelv2)>3
p2=polyfit(sinelv2,SNR2,2);
y2=polyval(p2,sinelv);
SNR2dt=SNR2a-y2;
if allsats(aa)<33
L2car=(299792458/1227.60e06); % GPS
else
L2car=(299792458/1246e06); % GLONASS
end
maxf2=numel(sinelv)/(2*(max(sinelv)-min(sinelv)));
ovs2=round(L2car/(2*prec1*(max(sinelv)-min(sinelv))));
tmpsnr2=SNR2dt;
tmpsnr2(isnan(tmpsnr2)==1)=[];
if size(tmpsnr2,1)>0
[psd2,f2]=plomb(SNR2dt,sinelv,maxf2,ovs2,'normalized');
reflh2=f2.*0.5*L2car;
[~,id2]=max(psd2(:));
pks2=findpeaks(psd2);
pks2=sort(pks2);
if reflh2(id2) < ahgt-ahgt_bounds || reflh2(id2) > ahgt+ahgt_bounds...
        || max(psd2)<10*mean(pks2(1:end-1))
    nol2=1;
end
else
    nol2=1;
end
else
    nol2=1;
end

if tempsnr==1 || templsp==1
width = 3.6;     % Width in inches % 3.5 was for putting two on one line i think
height = 1.3; % was 1.3    % Height in inches
fsz = 11;      % Fontsize
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
%left=0;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
close all
end
%%%%%%%%%%%%%%%
% SNR PLOT
if tempsnr==1 && nol1==0
figure('visible','off')
plot(asind(sinelv),SNR1dt,'b','linewidth',0.8)
xlabel('Elevation angle','interpreter','latex','fontsize',fsz)
ylabel('$\delta SNR$','interpreter','latex','fontsize',fsz)
set(gca,'ticklabelinterpreter','latex','fontsize',fsz-1)
axis([elv_low elv_high -200 200])
print([pwdstr,'/tempoutput/SNR_',num2str(allsats(aa)),'_',...
    num2str(round(mean(sat2tr(:,4)))),...
    '.png'],'-dpng','-r300')
end
% PERIODOGRAM
if templsp==1 && nol1==0
figure('visible','off')
plot(reflh1,psd,'r')
axis([0 inf 0 30])
xlabel('Frequency','interpreter','latex','fontsize',fsz)
ylabel('Power','interpreter','latex','fontsize',fsz)
set(gca,'ticklabelinterpreter','latex','fontsize',fsz-1)
print([pwdstr,'/tempoutput/LSP_',num2str(allsats(aa)),'_',...
    num2str(round(mean(sat2tr(:,4)))),...
    '.png'],'-dpng','-r300')
end

if tempfresnel==1 %&& nol1==0
    % currently for L1 signal, could change to L2, see function description
    fresout=[pwdstr,'/tempoutput'];
    googlefresnel(plat,plon,sat2tr(1,2),sat2tr(1,3),doo,fresout,reflh1(id),1)
    googlefresnel(plat,plon,sat2tr(end,2),sat2tr(end,3),doo,fresout,reflh1(id),1)
end

cnt=cnt+1;
if ltdrun~=0
if cnt==ltdrun
    return
end
end
slvlr(cnt,1)=tdatenum+mean(sat2tr(:,4))/86400;          % datenum
slvlr(cnt,2)=doo;                                       % Sat #
if tropd~=0                                            % tan(th)/dth/dt
slvlr(cnt,3)=tand(mean(thetarefr(:)))/(((pi/180)*(thetarefr(end)-thetarefr(1)))...
    /(sat2tr(end,4)-sat2tr(1,4)));
clear thetarefr
else
slvlr(cnt,3)=tand(mean(sat2tr(:,2)))/(((pi/180)*(sat2tr(end,2)-sat2tr(1,2)))...
    /(sat2tr(end,4)-sat2tr(1,4)));
end
slvlr(cnt,4)=min(sat2tr(:,2));                          % THETA MIN
slvlr(cnt,5)=max(sat2tr(:,2));                          % THETA MAX
slvlr(cnt,6)=nanmean(sat2tr(:,3));                      % MEAN AZI
% signal stats
% l1
if nol1==0
slvlr(cnt,7)=reflh1(id);                                % L1 slvl
else
slvlr(cnt,7)=NaN;
end
slvlr(cnt,8)=nanmean(y1);                               % mean mag. tSNR
slvlr(cnt,9)=max(psd);                                  % the peak
slvlr(cnt,10)=var(SNR1dt);                              % the variance
% l2
if nol2==0
slvlr(cnt,11)=reflh2(id2);                              % L2 slvl
slvlr(cnt,12)=nanmean(y2);                              % mean mag. tSNR
slvlr(cnt,13)=max(psd2);                                % the peak
slvlr(cnt,14)=var(SNR2dt);                              % the variance
else
slvlr(cnt,11:14)=NaN;
end 

lspy_flim=300;
lspy(cnt,:)=interp1(f,psd,1:1:lspy_flim);


end
end

if cnt==0
    slvlr=[];
    lspy=[];
end

end

