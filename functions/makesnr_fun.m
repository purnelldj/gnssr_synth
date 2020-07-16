function [snr_data,slvlr,lspy] = makesnr_fun(station,tdatenum,slvlrobs,tgstring,sp3str,effects,tempsnr,templsp)

%%

% This function creates synthetic SNR data as per reference:
% Purnell et al. (2020) "Quantifying the Uncertainty in Ground-Based 
% GNSS-Reflectometry Sea Level Measurements" (IEEE JSTARS)
%
% INPUTS
%
% station: 4 letter character vector for identifier of GNSS station e.g., 'sc02'
% tdatenum: day in matlab datenum format
% slvlrobs: observed stats from observed SNR data, output from
% tgstring: character vector with path to file of sea level variations
% should be in format of xaxis (datenum array), slvl (sea level array in m)
% sp3string: character vector with path to sp3 orbit file
% effects: (1,4) double array to turn effects on (1) or off (0)
% the first entry dynamic surface
% the second entry is noise
% the third entry is tropospheric delay
% the fourth entry is surface roughness
% tempsnr: choose 1 to output figures of detrended synthetic SNR data to the
% directory 'data/output_figs'
% templsp: choose 1 to output figures of Lomb-Scargle Periodogram data to the
% directory 'data/output_figs'
% ltdrun: use this as an option to produce a set amount of output figures
% and then stop
%
% OUTPUTS
% 
% snr_data: 7 column output of synthetic SNR data
% in the same format as output from RinexSNRv2
% column 1 is satellite PRN
% column 2 is elevation angle (degrees)
% column 3 is azimuth angle (degrees)
% column 4 is time in seconds of day
% column 7 is the synthetic SNR data in dB-Hz (main output)
% slvlr: a (k,18) array of stats of synthetic SNR data in same format as
% slvlr output from 'analyzesnr_fun.m' (see function description)
% lspy: a (k,300) array of LSP spectrums for k satellite arcs
% the frequency is interpolated by default from 1:300, this can be changed
% by setting the value of lspy_flim

ltdrun=0; % change if you just want to produce a set number of figures

curdt=datetime(tdatenum,'convertfrom','datenum');
curjd=juliandate(curdt);

% it might be important to clear these variables if you're calling the
% function in a loop / doing multiple days
clear slvlr lspy lspyt snr_data

% for doing multiple runs simultaneously, making sure that noise is not the same
%warning('off','all')
rng('shuffle')

% these should be your local directory (full strings)
pwdstr=pwd;
functionstr=[pwdstr,'/functions']; 
% directory that contains GTP2w, mpsim, lla2ecef,
% station_inputs
% there is also a matlab function lla2ecef that is probably better
realsignal=1; % turn on (1) or off (0) using stats from real signal

dyncor=effects(1);
addnoise=effects(2);
tropdelay=effects(3);
sfr=effects(4);

% other
nonstat_sfc=1;
glonass=0; % add glonass satellites
reslnfix=0; % SNR resolution in dB-Hz
mag_adj=1; % turn on (1) or off (0)
fixed_mag=0; % set the value (in dB-Hz, or 0 if taking from data
% should only be INSTEAD of tropdelay

% don't change these
numpks=15; % choose the number of peaks / oscillations to add for noise
doelvlims=0;
setphase=180; % change phase of SNR signal, increase to move to the right
maxf1fix=0; % to stop
hgtfix=0;
scaletide=1; % scale tides by a factor !!LEAVE AT 1!!
signal='L1'; % code needs to be modified to use other signals

addpath(functionstr)
addpath([functionstr,'/mpsim'])
run([functionstr,'/station_inputs/',station,'_input'])
delete([pwdstr,'/tempoutput/*'])

if realsignal==0
    resln=0;
end
if sfr==0
    seasfc_rough=0;
end
if reslnfix~=0
   resln=reslnfix;
end
if hgtfix~=0
    ahgt=hgtfix;
end

lat=plat;
lon=plon;

% NOTE: FIND ANTENNAS AVAILABLE BY LOOKING IN
% mpsim/lib/snr/data/ant/
% antenna_mod='TRM29659.00'; % choose antenna model
% radome_mod='SCIT'; % choose radome model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART ONE: GET ORBIT DATA

dt=abs(dt);

% tropd
if tropdelay==1
curjd2=curjd-2400000.5;
[~,~,hell]=ecef2lla(staxyz(1),staxyz(2),staxyz(3)); % HELL!!!!!
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

% refraction
if tropdelay==2
curjd2=curjd-2400000.5;
[~,~,hell]=ecef2lla(staxyz(1),staxyz(2),staxyz(3)); % HELL!!!!!
hgtlim(1) = hell-ahgt-ahgt_bounds;
hgtlim(2) = hell-ahgt+ahgt_bounds;
dlat(1) = plat*pi/180.d0;
dlon(1) = plon*pi/180.d0;
nstat = 1;
it = 0;
[plim(1),tlim(1),~,~,~,~,~,~,~] = gpt2_1w (curjd2,dlat,dlon,hgtlim(1),nstat,it);
[plim(2),tlim(2),~,~,~,~,~,~,~] = gpt2_1w (curjd2,dlat,dlon,hgtlim(2),nstat,it);
end

% ORBIT DATA
if glonass==1
    satsn=32+24;
else
    satsn=32;
end
satmatr={'G01';'G02';'G03';'G04';'G05';'G06';'G07';'G08';'G09';...
    'G10';'G11';'G12';'G13';'G14';'G15';'G16';'G17';'G18';'G19';...
    'G20';'G21';'G22';'G23';'G24';'G25';'G26';'G27';'G28';'G29';...
    'G30';'G31';'G32';...
    'R01';'R02';'R03';'R04';'R05';'R06';'R07';'R08';'R09';...
    'R10';'R11';'R12';'R13';'R14';'R15';'R16';'R17';'R18';'R19';...
    'R20';'R21';'R22';'R23';'R24'}; %obs would need to add in galileo
orbit_data_all=[];
for satind=1:satsn

satname=char(satmatr(satind));
fid=fopen(sp3str,'r');

i=1;
tline=fgets(fid);
while i<23
    i=i+1;
    tline=fgets(fid);
end
    clear xyz
        
tt=0;
while ~feof(fid)
    tt=tt+1;
    sattime(tt)=datenum(str2double(tline(4:7)),str2double(tline(9:10)),str2double(tline(12:13)),...
        str2double(tline(15:16)),str2double(tline(18:19)),str2double(tline(21:31)));
    tline=fgets(fid);
    while strcmp(tline(1),'*')==0 && ~feof(fid)
        if strcmp(tline(2:4),satname)==1
            xyz(tt,:)=[str2double(tline(5:18)) str2double(tline(19:32)) str2double(tline(33:46))];
        end
        tline=fgets(fid);
    end
end
if exist('xyz')==0
    continue
end
temp=sum(xyz,2);
tempdel=temp(:)==0;
xyz(tempdel,:)=[];
sattime(tempdel)=[];
t=(sattime-sattime(1)).*86400;
t=t.';
ts=0:dt:86400-dt;
ts=ts.';
a=0;
interp(:,1)=spline(t,xyz(:,1),ts);
interp(:,2)=spline(t,xyz(:,2),ts);
interp(:,3)=spline(t,xyz(:,3),ts);
xyz=xyz.*1000;
interp=interp.*1000;
if exist('orbit_data')==1
    clear orbit_data
end 
ind=0;
for tt=1:size(ts,1)
[azi,elv]=gnss2azelv(staxyz,interp(tt,:),lat,lon);
if (azi>azi_low && azi<azi_high) && (elv>elv_low && elv<elv_high)
    if exist('extraazilims')==1
        if azi>extraazilims(1) && azi<extraazilims(2)
            continue
        end
    end
ind=ind+1;
orbit_data(ind,1)=satind;
orbit_data(ind,2)=elv;
orbit_data(ind,3)=azi;
orbit_data(ind,4)=ts(tt);
orbit_data(ind,5)=NaN; %SNR L1
end
end
clear interp

fclose(fid);

if ind>0
orbit_data_all(size(orbit_data_all,1)+1:size(orbit_data_all,1)+size(orbit_data,1),:)=...
    orbit_data(1:end,:);
clear orbit_data
end

end

orbit_data=orbit_data_all;
allsats=unique(orbit_data(:,1));

% NOW GET TIDE GAUGE DATA
load(tgstring)
datevecee=tdatenum+1;
checknans=isnan(slvl(:))==1;
if sum(checknans)>0
    disp('nans in sea level data')
    return
end
tidey=scaletide*slvl;
tttmp=tdatenum:dt/86400:datevecee-dt/86400;
tidedt=interp1(xaxis,tidey,tttmp,'linear');
dell=find(isnan(tidedt(:)));
if numel(dell)==1 && dell~=size(tidedt,2)
    tidedt(dell)=tidedt(dell+1);
end

clear xaxis slvl

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART TWO: GENERATE SNR DATA

snr_data=[];
cnt_st=0;
cnt=0;

for aa=1:numel(allsats)

repgaps=1;    
ijk=0;
while ijk~=repgaps
    
    dt=abs(dt);    
    ijk=ijk+1;
    if exist('sat2tr')==1
        clear sat2tr
    end
    
doo=allsats(aa);
tmpid=orbit_data(:,1)==doo;
sat2=orbit_data(tmpid,:);
in=sat2(:,3)>azi_low & sat2(:,3)<azi_high;
sat2trr=sat2(in,:);
in=sat2trr(:,2)<elv_high & sat2trr(:,2)>elv_low;
sat2tr=sat2trr(in,:);
if exist('azi_mask')==1
    out=sat2tr(:,3)>azi_mask(1) & sat2tr(:,3)<azi_mask(2);
    sat2tr(out)=[];
end
%%%%%%%
if hgtfix~=0
    ahgt=hgtfix;
end
times=sat2tr(:,4);
if size(times,1)<1
    continue;
end

% to sort satellites which have more than one valid overpass
gaps=diff(times(:,1))>dt*2;
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
if size(sat2tr,1)<2
    continue
end

fwd=0;
if sat2tr(2,2)-sat2tr(1,2)>0
    fwd=1;
end
for a=3:size(sat2tr,1) % so this is basically just deleting all data if it ...
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


if doelvlims==1
if abs(min(sat2tr(:,2))-elv_low)>elvlims || abs(max(sat2tr(:,2))-elv_high)>elvlims
    disp('elvlims')
    continue
end
end

% this is a minimum time limit!!!
if abs(sat2tr(end,4)-sat2tr(1,4))<300
    continue
end

if realsignal==1
    % now see if this was a recorded sea level measurement
    %return
    in=slvlrobs(:,2)==allsats(aa) & slvlrobs(:,1)>tdatenum+min(sat2tr(:,4))/86400 ...
    & slvlrobs(:,1)<tdatenum+max(sat2tr(:,4))/86400;
    slvlrt=slvlrobs(in,:);
    if size(slvlrt,1)<1
        disp('data not found')
        continue
    end
    if size(sat2tr,1)<300/dt
        disp('weird extra limit')
        continue
    end
    inelvs=sat2tr(:,2)>=slvlrt(4) & sat2tr(:,2)<=slvlrt(5);
    sat2tr=sat2tr(inelvs,:);
    % this is a minimum time limit!!!
    if size(sat2tr,1)<1
        disp('skip')
        continue
    end
end

% MOVING H DATA FROM TIDE GAUGE
if nonstat_sfc==1
tidedtt=interp1(0:dt:86400-dt,tidedt,sat2tr(:,4),'nearest'); % here going forwards in time (i hope)
h_temp=ahgt-tidedtt;
h_temp=h_temp.';
else
h_temp=ahgt*ones(size(sat2tr,1));
end

sinelv=sind(sat2tr(:,2)).';
theta=sat2tr(:,2).';

if tropdelay==2
hsfc=hell-h_temp;
psfc=interp1(hgtlim,plim,hsfc,'linear');
tsfc=interp1(hgtlim,tlim,hsfc,'linear');
dele = (1/60)*510*psfc/((9/5*tsfc+492)*1010.16).*cotd(theta+7.31./(theta+4.4));
theta=theta+dele;
sinelvd=sinelv;
sinelv=sind(theta);
sinelvrefr=sinelv;
thetarefr=theta;
end

if tropdelay==1
hsfc=hell-h_temp;
psfc=interp1(hgtlim,plim,hsfc,'linear');
tmsfc=interp1(hgtlim,tmlim,hsfc,'linear');
esfc=interp1(hgtlim,elim,hsfc,'linear');
clear thetarefr
for jj=1:numel(theta)
tau(jj)=trop_delay_tp(curjd2,plat,hell,hsfc(jj),theta(jj),pant,tmant,eant,ah,aw,lambda,psfc(jj),...
    tmsfc(jj),esfc(jj));
thetarefr(jj)=asind(sind(theta(jj))+0.5*tau(jj)/h_temp(jj));
end
theta=thetarefr;
sinelvd=sinelv;
sinelv=sind(thetarefr);
sinelvrefr=sinelv;
end

init
sett=snr_settings();

% SIGNAL
if allsats(aa)<33
sett.opt.gnss_name = 'gps';
L1car=(physconst('LightSpeed')/1575.42e06); % for GPS
else
sett.opt.gnss_name = 'glonass';
L1car=(299792458/1602e06); % for GLONASS
end
sett.opt.freq_name = signal;
% ANTENNA / INSTRUMENTATION
sett.ant.model = antenna_mod;
sett.ant.radome = radome_mod;
sett.ref.ignore_vec_apc_arp = true; % gets rid of offset between antenna 
%phase centre and reference point??

% PHASE
sett.bias.phase_interf = setphase;

% IMPORTANT
if fwd==1
input_hgt = h_temp(1);
else
input_hgt = h_temp(end);
end

if dyncor==0
    input_hgt=mean(h_temp);
end

% NON STATIC SEA SURFACE
if dyncor==1
h=h_temp;
if fwd==0
    % at this point from high to low elv
    h=fliplr(h);
    sinelv=fliplr(sinelv);
    theta=fliplr(theta);
end
clear dsinelv_ns dhdt_ns dthdt_ns dsinelv_s theta_s
theta_s(1)=theta(1);
for i=2:numel(h)
    dsinelv_ns(i)=sinelv(i)-sinelv(i-1);
    dhdt_ns(i)=(h(i)-h(i-1))/dt;
    dthdt_ns(i)=(pi/180)*(theta(i)-theta(i-1))/dt;
    dsinelv_s(i)=dsinelv_ns(i)*(mean(h(i-1:i))+dhdt_ns(i)*...
        tand(mean(theta(i-1:i)))/dthdt_ns(i))/input_hgt;
    theta_s(i)=asind(sind(theta_s(i-1))+dsinelv_s(i));
end
thetass=theta_s(1);
thetase=theta_s(end);
sinelv_stat=sind(theta_s);
else
    thetass=min(theta);
    thetase=max(theta);
end

sett.ref.height_ant = input_hgt;

sett.sat.elev_lim = [thetass,thetase];
% make theta evenly spaced, which is fine assuming dth/dt is fixed
sett.sat.num_obs = 1000; 

% REFLECTING SURFACE
sett.sfc.bottom_material = 'seawater';
% NOTE BELOW OPTIONS DIDN'T CHANGE ANYTHING
%temp = 8; % extremes probably 8 to 15 ?? 2 is default
%salinity = 30; % extremes are probably 24 to 28....26 should be default
%perm=get_permittivity(struct('name','seawater variable','temperature',...
    %num2str(temp),'salinity',num2str(salinity)));

% SEA SURFACE ROUGHNESS
sett.sfc.height_std = seasfc_rough; 

setup = snr_setup (sett);
result = snr_fwd (setup);

% INTERPOLATING HERE TO MATCH ORBIT
snr_preproc=result.snr_db;
snr_preproc=snr_preproc.';
theta_init=result.sat.elev;

if dyncor==1
snr_preproc=interp1(theta_init,snr_preproc,theta_s);
else
    if fwd==1
snr_preproc=interp1(theta_init,snr_preproc,theta);
    else
snr_preproc=interp1(theta_init,snr_preproc,fliplr(theta));        
    end
end

dell=isnan(snr_preproc)==1;
snr_preproc(dell)=[];
if fwd==0
    sat2tr=flipud(sat2tr);
end
sat2tr(dell,:)=[];
if tropdelay==2
sinelvrefr(dell)=[];
end
sinelv=sind(sat2tr(:,2)).';
theta=sat2tr(:,2).';
if fwd==0
    sat2tr=flipud(sat2tr);
end

% MEAN MAGNITUDE OF SIGNAL
if mag_adj==1
snr_preproc_beformag=snr_preproc;
SNRt=sqrt(10.^(snr_preproc./10));
pt=polyfit(sinelv,SNRt,2);
yt=polyval(pt,sinelv);
if fixed_mag==0
diffmag=nanmean(yt)-slvlrt(8);
else
diffmag=nanmean(yt)-fixed_mag;  
end
SNRtt=SNRt-diffmag;
snr_preproc=10.*log10(SNRtt.^2);
end

% NOW ADD NOISE
if addnoise==1
    redo=1;
    iter=0;
    SNRt=sqrt(10.^(snr_preproc./10));
    pt=polyfit(sinelv,SNRt,2);
    yt=polyval(pt,sinelv);
    SNRdt=SNRt-yt;
    maxf1=numel(sinelv)/(2*(max(sinelv)-min(sinelv)));
    if maxf1fix~=0
        maxf1=maxf1fix;
    end
    prec1=0.001;
    ovs=round(L1car/(2*prec1*(max(sinelv)-min(sinelv))));
    ovs=round(ovs/10);
    %[psd,f]=plomb(SNRdt,sinelv,maxf1,round(ovs/10),'normalized'); %
    fi=1:1:maxf1;
    [psd,f,~,~,~,~]=fLSPw(sinelv,SNRdt,fi,0.05,ovs);
    psd=cell2mat(psd);
    f=cell2mat(f);
    [~,id]=max(psd);
    prenoisepkf=f(id);
    psdt=2*sqrt(psd.*2*var(SNRdt)/numel(sinelv));
    pkinit=max(psdt);
    killow=5;
    while redo==1
    iter=iter+1;
    freqsin=killow+rand(numpks,1)*(maxf1-killow);
    if strcmp(station,'pbay')==1
        freqsin=killow+rand(numpks,1)*(300-killow); % % NOT TO MAX(f)
    end
    phaseoff=rand(numpks,1)*2*pi-pi;
    pkin=slvlrt(9);
    varint=slvlrt(10);
    powerin=rand(numpks,1);
    coefsin=[0.1;pkin/pkinit*10]; % was 0.001
    [powerout]=scalenoisepower(station,powerin,freqsin);
    coefsin=log(coefsin);
    tempfun=@(coefs) noisefunresid(coefs,snr_preproc,powerout,freqsin,sinelv,phaseoff,resln,...
        pkin,varint,maxf1,ovs);
    options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','off');
    coefs_ls=lsqnonlin(tempfun,coefsin,[],[],options); %lsqnonlin or fsolve??
    coefs=coefs_ls;
    [snr_recr,yt2]=noisefun(coefs,snr_preproc,powerout,freqsin,sinelv,phaseoff,resln);
    snr_recr=real(snr_recr); % this happened at sc02 for modt doy 220
    yt2=real(yt2);
    fi=1:1:maxf1;
    [psdn,fn,A2,~,~,~]=fLSPw(sinelv,snr_recr,fi,0.05,ovs);
    psdn=cell2mat(psdn);
    fn=cell2mat(fn);
    varresid=var(snr_recr)-varint;
    pkresid=max(psdn)-pkin;
    residd=[varresid;pkresid];
    reflh1=fn.*0.5*L1car;
    [~,id]=max(psdn);
    pks=findpeaks(psdn);
    pks=sort(pks);
    if ~isnan(slvlrt(7))
    if reflh1(id) > ahgt-ahgt_bounds && reflh1(id) < ahgt+ahgt_bounds...
        && max(psdn)>10*mean(pks(1:end-1)) && sum(abs(residd))<1e-3 && abs(fn(id)-prenoisepkf)<5
        redo=0;
    end
    elseif reflh1(id) < ahgt-ahgt_bounds || reflh1(id) > ahgt+ahgt_bounds...
        || max(psdn)<10*mean(pks(1:end-1))
        redo=0;
    end
    if iter>=1 && redo==1
        disp('iterating')
        if iter>20
        disp('too many times - skipped')
        redo=0;
        end
    end
    end
    pkout=max(psdn);
    varout=var(snr_recr);
    snr_recr=snr_recr+yt2;
    snr_preproc=10.*log10(snr_recr.^2);
end
clear residd

% NOW IMPOSE SNR RESOLUTION
if resln~=0
snr_init=snr_preproc;
mods=mod(snr_init,resln);
for i=1:size(mods,2)
    if mods(i)>=resln/2
        snr_resln(i)=snr_init(i)-mods(i)+resln;
    else
        snr_resln(i)=snr_init(i)-mods(i);
    end
end
snr_preproc=snr_resln;
clear snr_resln
end

SNR1=sqrt(10.^(snr_preproc./10));
p1=polyfit(sinelv,SNR1,2);
y1=polyval(p1,sinelv);
SNR1dt=SNR1-y1;
%SNR1dt=SNR1-yt;
maxf1=numel(sinelv)/(2*(max(sinelv)-min(sinelv)));
prec1=0.001;
ovs=round(L1car/(2*prec1*(max(sinelv)-min(sinelv))));
%[psd,f]=plomb(SNR1dt,sinelv,maxf1,ovs,'normalized'); %
fi=1:1:maxf1;
[psd,f,A2,~,~,~]=fLSPw(sinelv,SNR1dt,fi,0.05,ovs);
psd=cell2mat(psd);
f=cell2mat(f);
reflh1=f.*0.5*L1car;
[~,id]=max(psd(:));

testresids=0;
if testresids==1 && dyncor==1
    corrr=mean(dhdt_ns(1:end))*tand(mean(sat2tr(:,2)))/...
        (((pi/180)*(sat2tr(end,2)-sat2tr(1,2)))/(sat2tr(end,4)-sat2tr(1,4)));
    if fwd==0
        corrr=-corrr; % because of the dhdt
    end
[psds,fs,A2,~,~,~]=fLSPw(sind(theta_s),SNR1dt,fi,0.05,ovs);
psds=cell2mat(psds);
fs=cell2mat(fs);
reflh1s=fs.*0.5*L1car;
[~,ids]=max(psds(:));
resid_s=input_hgt-reflh1s(ids);
residcor=abs(mean(h_temp)-(reflh1(id)-corrr)-resid_s);
residnocor=abs(mean(h_temp)-reflh1(id)-resid_s);
improvement=residnocor-residcor;
disp(['improvement=',num2str(improvement)])
if improvement<-0.05
    return
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
snrsave=10.*log10(SNR1.^2);
if fwd==0
    snrsave=fliplr(snrsave);
end
sat2tr(:,7)=snrsave;
sat2tr(:,5)=NaN;
sat2tr(:,6)=NaN;
sat2tr(:,8)=NaN;
snr_data=[snr_data;sat2tr];
%%%%%%%%%%%%%%%%%%%%%%%%%%%

nol1=0;
nol2=1;
pks=findpeaks(psd);
pks=sort(pks);
if reflh1(id) < ahgt-ahgt_bounds || reflh1(id) > ahgt+ahgt_bounds...
        || max(psd)<10*mean(pks(1:end-1)) || isnan(slvlrt(7)) %|| sum(abs(resid))>1e-3
    nol1=1;
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
if tempsnr==1
figure('visible','off')
plot(sinelv,SNR1dt,'b','linewidth',0.8)
xlabel('$\sin{\theta}$','interpreter','latex','fontsize',fsz)
ylabel('$\delta SNR$','interpreter','latex','fontsize',fsz)
set(gca,'ticklabelinterpreter','latex','fontsize',fsz-1)
axis([sind(elv_low) sind(elv_high) -150 150])
print([pwdstr,'/tempoutput/SNR_',num2str(allsats(aa)),'_',...
    num2str(round(mean(sat2tr(:,4)))),...
    '.png'],'-dpng','-r300')
end
% PERIODOGRAM
if templsp==1
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

cnt=cnt+1;
if ltdrun~=0
if cnt==ltdrun
    return
end
end
slvlr(cnt,1)=tdatenum+mean(sat2tr(:,4))/86400;          % datenum
slvlr(cnt,2)=doo;                                       % Sat #
if tropdelay~=0                                            % tan(th)/dth/dt
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

end

