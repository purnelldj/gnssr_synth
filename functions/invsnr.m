function [sfacsjs,sfacspre,hinit,xinit,consts_out,roughness] = invsnr(tdatenum,station,snrdir,kspac,tlen,decimate,...
    satconsts,sfcrough,altelvlims,largetides)

% this code is for inverse modelling of SNR data to get GNSS-R measurements
% as per Strandberg et al. (2016)

% INPUTS
% tdatenum: day or time in datenum format
% station: station identifier string
% snrdir: path to SNR data, should also be in cellstr format
% kspac: average node spacing in days (e.g., 2/24 is 2 hours)
% tlen: length of window for analysis in days, split into 3 and the middle
% period is saved (e.g., setting to 3 means that the middle day is saved)
% decimate: decimate input SNR data (in seconds)
% satconsts: 1 by 3 double, 1 or 0 to include satellite constellations
% GPS, GLONASSS, GALILEO (e.g., [1 0 1] for GPS and GALILEO)
% sfcrough: initial sfc roughness guess ('s') in metres for least squares
% adjustment
% altelvlims: to overwrite the elevation limits in the station input file
% and use alternate ones, e.g., [5 15] for 5 to 15 degrees
% largetides: changes the initial guess for node values (1 or 0)

% OUTPUTS
% these outputs are to go with the function 'invsnr_plot.m'
% sfacsjs: node values estimated from Strandberg et al. analysis
% sfacspre: node values estimated by using spectral analysis adjustment
% hinit: initial spectral analysis estimates
% xinit: timing of spectral analysis estimates
% consts_out: the other variables estimated as part of least squares
% adjustment
% roughness: the roughness parameter that is estimated in least squares
% adjustment

p=2; % bspline order
stdfac=3;

gps=satconsts(1);
glo=satconsts(2);
gal=satconsts(3);

pwdstr=pwd;
addpath([pwdstr,'/functions/'])
addpath([pwdstr,'/functions/bspline'])
run([pwdstr,'/functions/station_inputs/',station,'_input'])
if glo==1
    load('glonasswlen.mat')
end
if decimate~=0
    dt=decimate;
end
if numel(altelvlims)>0
    elv_low=altelvlims(1);
    elv_high=altelvlims(2);
end

tdatenum=tdatenum-tlen/3;
 
curdt=datetime(tdatenum+tlen/3,'convertfrom','datenum');
disp(char(curdt))

% work out if need two days or just one
mlen=1;
if tdatenum+tlen-mod(tdatenum+tlen,1)-(tdatenum-mod(tdatenum,1))>0
    mlen=tdatenum+tlen-mod(tdatenum+tlen,1)-(tdatenum-mod(tdatenum,1))+1;
    if mod(tdatenum+tlen,1)==0
        mlen=mlen-1;
    end
end

% AT THIS POINT THE tdatenum AND CURJD ARE THE SAME

% try getting day in order of station first, then after detrended, organise
% by time, then that's it
snrfile=[];
for ll=1:numel(snrdir)
    snrfilet=[];
for m=1:mlen
    % NOW LOADING FROM tdatenum ONWARDS
    tdatenumt=tdatenum+m-1;
    curdtt=datetime(tdatenumt,'convertfrom','datenum');
    strdayt=char(datetime(curdtt,'format','DDD'));
    stryrst=char(datetime(curdtt,'format','yy'));
    if exist([char(snrdir(ll)),'/',num2str(tdatenumt),'.mat'],'file')==2
        load([char(snrdir(ll)),'/',num2str(tdatenumt),'.mat'])
    elseif exist([char(snrdir(ll)),'/',station,strdayt,'0.',stryrst,'snr'])==2
        snr_data=dlmread([char(snrdir(ll)),'/',station,strdayt,'0.',stryrst,'snr']);
    else
        disp('missing data')
        miss=1;
        break
    end
    snr_data(:,9)=tdatenum-mod(tdatenum,1)+m-1+snr_data(:,4)./86400;
    snr_data(:,10)=ll;
    % now get rid of data outside window
    tmpout=snr_data(:,9)<tdatenum | snr_data(:,9)>=tdatenum+tlen;
    snr_data(tmpout,:)=[];
    snrfilet=[snrfilet;snr_data];
end
if exist('miss')~=0
    break
end
snrfile=[snrfile;sortrows(snrfilet,1)];
clear snrfilet snr_data
end

if decimate~=0
    modspl=mod(snrfile(:,4),decimate);
    delt=modspl(:)>0;
    snrfile(delt,:)=[];
end

in=snrfile(:,3)>azi_low & snrfile(:,3)<azi_high;
snrfile=snrfile(in,:);
in=snrfile(:,2)<elv_high & snrfile(:,2)>elv_low;
snrfile=snrfile(in,:);
if exist('azi_mask')==1
    out=snrfile(:,3)>azi_mask(1) & snrfile(:,3)<azi_mask(2);
    snrfile(out,:)=[];
end
tmp=snrfile(:,7)==0;
snrfile(tmp,7)=NaN;
tmp=snrfile(:,8)==0;
snrfile(tmp,8)=NaN;
snrfile(:,7)=sqrt(10.^(snrfile(:,7)./10));
snrfile(:,8)=sqrt(10.^(snrfile(:,8)./10));
dtdv=dt/86400;

% getting desired satellite constellations
if gps==0
    delete=snrfile(:,1)<32+1;
    snrfile(delete,:)=[];
end
if glo==0
    delete=snrfile(:,1)>32 & snrfile(:,1)<32+24+1;
    snrfile(delete,:)=[];
end
if gal==0
    delete=snrfile(:,1)>32+24;
    snrfile(delete,:)=[];
end

%%% NEED TO DETREND SNR DATA
ind=1;
cursat=snrfile(1,1);
stind=1;
stopp=1;
s1ind=0;
sinelv1_all=[];
snr1_all=[];
t1_all=[];
satno_all=[];
antno_all=[];
prec1=0.001;
if snrfile(2,2)-snrfile(1,2)<0
    fwd2=0;
else
    fwd2=1;
end
while stopp==1
    ind=ind+1;
    if ind==size(snrfile,1)
        stopp=0;
    end
    if snrfile(ind,2)-snrfile(ind-1,2)<0
        fwd1=0;
    else
        fwd1=1;
    end
    curdt=snrfile(ind,9)-snrfile(ind-1,9);
    if ind-stind>1
    if snrfile(ind,1)~=cursat || ind==size(snrfile,1) || abs(curdt)>=3*dtdv...
            || fwd2~=fwd1 || snrfile(ind,2)==snrfile(ind-1,2)
        if ind-stopp-stind<(300/dt)
            snrfile(stind:ind-stopp,:)=[];
            ind=stind;
        else
        sinelvt1=snrfile(stind:ind-stopp,2);
        snr1tmp=snrfile(stind:ind-stopp,7);
        t1tmp=snrfile(stind:ind-stopp,9);
        satnotmp=snrfile(stind:ind-stopp,1);
        antnotmp=snrfile(stind:ind-stopp,10);
        del=isnan(snr1tmp(:,1))==1;
        sinelvt1(del,:)=[];
        if numel(sinelvt1)>2
        sinelvt1=sind(sinelvt1);
        snr1tmp(del,:)=[];
        t1tmp(del,:)=[];
        satnotmp(del,:)=[];
        antnotmp(del,:)=[];
        p1=polyfit(sinelvt1,snr1tmp,2);
        y1=polyval(p1,sinelvt1);
        sinelv1_all=[sinelv1_all;sinelvt1];
        snr1tmp=snr1tmp-y1;
        snr1_all=[snr1_all;snr1tmp];
        t1_all=[t1_all;t1tmp];
        satno_all=[satno_all;satnotmp];
        antno_all=[antno_all;antnotmp];
        % TO GET INITIAL H TIME SERIES
        if sinelvt1(2,1)-sinelvt1(1,1)<0
            sinelvt1=flipud(sinelvt1);
            snr1tmp=flipud(snr1tmp);
        end
        if snrfile(stind,1)<33
        L1car=(299792458/(1575.42e06/1.023)); % for GPS
        elseif snrfile(stind,1)<57
        L1car=glonasswlen(snrfile(stind,1)-32);
        elseif  snrfile(stind,1)>56
        L1car=(299792458/(1575.42e06/1.023));
        end
        maxf1=numel(sinelvt1)/(2*(max(sinelvt1)-min(sinelvt1)));
        ovs=round(L1car/(2*prec1*(max(sinelvt1)-min(sinelvt1))));
        if sum(diff(sinelvt1)>0)==numel(sinelvt1)-1 || sum(diff(sinelvt1)>0)==0
        [psd,f]=plomb(snr1tmp,sinelvt1,maxf1,ovs,'normalized'); %
        reflh1=f.*0.5*L1car;
        [~,id]=max(psd(:));
        pks=findpeaks(psd);
        pks=sort(pks);
        skiphere=0;
        else
        skiphere=1;
        end
        if skiphere==0
        if reflh1(id) > ahgt-ahgt_bounds && reflh1(id) < ahgt+ahgt_bounds
        s1ind=s1ind+1;
        hinit(s1ind)=reflh1(id);
        xinit(s1ind)=mean(t1tmp);
        snr_datatmp=snrfile(stind:ind-stopp,:);
        tanthter(s1ind)=tand(mean(snr_datatmp(:,2)))/(((pi/180)*...
            (snr_datatmp(end,2)-snr_datatmp(1,2)))...
            /((snr_datatmp(end,9)-snr_datatmp(1,9))*86400));
        siteinit(s1ind)=snr_datatmp(end,10);
        end
        end
        end
        end
        stind=ind;
        if ind<size(snrfile,1)
            cursat=snrfile(ind,1);
        end
    end
    end
    fwd2=fwd1;
end

indt=t1_all(:)>tdatenum+tlen/3 & t1_all(:)<tdatenum+2*tlen/3;
t1_allt=t1_all(indt);
maxt1gap=max(diff(sort(t1_allt)));
if maxt1gap>kspac
    disp('gap in data bigger than node spacing')
    sfacsjs=NaN;
    sfacspre=NaN;
    hinit=NaN;
    xinit=NaN;
    consts_out=NaN;
    roughness=NaN;
    return
end

% for adjusting heights
if numel(snrdir)>1
    for ll=1:numel(snrdir)
        in=siteinit(:)==ll;
        meanhgts(ll)=nanmean(hinit(in));
        if ll>1
            hinit(in)=hinit(in)+meanhgts(1)-meanhgts(ll);
        end
    end
else
    meanhgts=0;
end

tmpinit=[xinit.' hinit.' tanthter.' siteinit.'];
tmpinit=sortrows(tmpinit,1);
xinit=tmpinit(:,1);
hinit=tmpinit(:,2);
tanthter=tmpinit(:,3);
siteinit=tmpinit(:,4);
hsmooth=smoothdata(hinit,'movmean',5);
diff1=abs(hsmooth-hinit);
std1=std(diff1);
delete=diff1(:,1)>stdfac*std1; % 2 for 4 stations
hinit(delete)=[];
xinit(delete)=[];
tanthter(delete)=[];
siteinit(delete)=[];

knots=[tdatenum*ones(1,p) ...
    tdatenum:kspac:tdatenum+tlen ...
    (tdatenum+tlen)*ones(1,p)];
nsfac=tlen/kspac+p;
sfacs_0=ahgt*ones(1,nsfac);
tempfun_init=@(sfacs) bspline_spectral(sfacs,p,knots,tanthter,xinit,1)-hinit.';
options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
    'Display','off'); % off for now but maybe should check sometimes
sfacs_init=lsqnonlin(tempfun_init,sfacs_0,[],[],options); %lsqnonlin or fsolve??
in=xinit(:)>tdatenum+tlen/3 & xinit(:)<tdatenum+2*tlen/3;
xinit=xinit(in).';
hinit=hinit(in).';
sfacspre=sfacs_init;

disp('joakim part')
doprev=1;
if exist('sfacsjs')==0 || doprev==0
if largetides==1
sfacs_0=sfacs_init;
else
sfacs_0=median(hinit)*ones(size(sfacs_init));
end
else
inds=tlen/(3*kspac)+2;
sfacs_0=sfacsjs(inds:end-1);
sfacs_0=[sfacs_0(1)*ones(1,p-1) sfacs_0 sfacs_0(end)*ones(1,tlen/(3*kspac)+p-1)];
end
consts=gps+glo+gal;
sfacs_0=[sfacs_0 zeros(1,consts*2)];
sfacs_0=[sfacs_0 sfcrough]; % OK SO SFR ROUGH SHOULD BE 0 IF MODEL OR 0.001
tempfun=@(sfacs) bspline_js(sfacs,t1_all,sinelv1_all,snr1_all,knots,...
    p,satno_all,gps,glo,gal,antno_all,meanhgts);
options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
    'Display','off');
tic
sfacs_ls=lsqnonlin(tempfun,sfacs_0,[],[],options); %lsqnonlin or fsolve??
toc
disp('****least squares done****')

% to plot model vs real snr data
plotmodsnr=0;
if plotmodsnr==1
    disp('plotting mod snr output')
residout=bsp_snrout(sfacs_ls,t1_all,sinelv1_all,snr1_all,knots,...
    p,satno_all,gps,glo,gal,antno_all,meanhgts,dtdv,elv_low,elv_high);
disp('done')
end

sfacsjs=sfacs_ls(1:end-consts*2-1);
consts_out=sfacs_ls(end-consts*2:end-1);
roughness=sfacs_ls(end);


end

