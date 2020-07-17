function [t_rh,rh_invjs,rh_invpre,rms_js,rms_pre] = invsnr_plot(startdate,enddate,invdir,kspac,tlen,...
    plotl,tgstring,makefig,roughnessplot)

%%

% this codes takes the output from invsnr.m and plots it and compares with
% a tide gauge (if there is one)

% INPUTS
% startdate: in datenum format
% enddate: in datenum format
% invdir: directory that contains invsnr.m output to plot
% kspac: average node spacing in days (e.g., 2/24 is 2 hours)
% tlen: length of window for analysis in days, split into 3 and the middle
% period is saved (e.g., setting to 3 means that the middle day is saved)
% plotl: the frequency of output time series, in days (e.g., 1/24 is
% hourly)
% tgstring: string of path to tide gauge data (or [] if no tide gauge data)
% makefig: set to 1 if you want to plot / save a figure
% roughnessplot:m add another panel with a daily plot of the roughness
% parameter

% OUPUTS
% t_rh: time vector in datenum format
% rh_invjs: output of b-spline reflector heights from Strandberg et al. analysis
% rh_invpre: output of b-spline reflector heights using spectral analysis
% rms_js: rms between tide gauge and invjs output
% rms_pre: rms between tide gauge and invdp output

pwdstr=pwd;
addpath([pwdstr,'/functions/bspline'])

meanormedian=2;
domeanheights=1;
p=2;

knots=[(startdate-1)*ones(1,p) ...
    startdate-1:kspac:enddate+2 ...
    (enddate+2)*ones(1,p)];
t_rh=startdate:plotl:enddate+1;

sfacspreall=[];
sfacsjsall=[];
x_init=[];
h_init=[];
roughness_all=[];
for jj=1:numel(invdir)
tdatenum=startdate-1;
while tdatenum<enddate
    tdatenum=tdatenum+1;
    inds=tlen/(3*kspac)+2;
    inde=(2*tlen)/(3*kspac)+1;
    if tdatenum==startdate
        inds=1;
    end
    if tdatenum==enddate
        inde=tlen/kspac+2*p/2;
    end
    if exist([char(invdir(jj)),'/',num2str(tdatenum),'.mat'],'file')==2
    load([char(invdir(jj)),'/',num2str(tdatenum),'.mat'])
    sfacspreall=[sfacspreall sfacspre(inds:inde)];
    sfacsjsall=[sfacsjsall sfacsjs(inds:inde)];
    x_init=[x_init xinit];
    h_init=[h_init hinit];
    roughness_all=[roughness_all roughness];
    else
        disp('missingdata')
    sfacspreall=[sfacspreall NaN(1,inde-inds+1)];
    sfacsjsall=[sfacsjsall NaN(1,inde-inds+1)];
    roughness_all=[roughness_all NaN];
    end
end    
end

if numel(invdir)>1
    indfac=numel(sfacsjsall)/numel(invdir);
    if domeanheights==1
    for kk=1:numel(invdir)
        % for adjusting heights
        meanhgts(kk)=nanmean(sfacsjsall((kk-1)*indfac+1:kk*indfac));
        if kk>1
            sfacsjsall((kk-1)*indfac+1:kk*indfac)=...
                sfacsjsall((kk-1)*indfac+1:kk*indfac)+meanhgts(1)-meanhgts(kk);
        end
    end
    end
    %meanhgts-meanhgts(1)
    %return
    for ii=1:indfac
        jj=0;
        dptmp=[];
        jstmp=[];
        while jj<numel(invdir)
            jj=jj+1;
            dptmp=[dptmp sfacspreall(ii+(jj-1)*indfac)];
            jstmp=[jstmp sfacsjsall(ii+(jj-1)*indfac)];
        end
        if meanormedian==1
        sfacsprep(ii)=nanmean(dptmp);
        sfacsjsp(ii)=nanmean(jstmp);
        else
        sfacsprep(ii)=nanmedian(dptmp);
        sfacsjsp(ii)=nanmedian(jstmp);
        end
    end
    roughlen=enddate-startdate;
    for ii=1:roughlen
        jj=0;
        roughtmp=[];
        while jj<numel(invdir)
            jj=jj+1;
            roughtmp=[roughtmp roughness_all(ii+(jj-1)*roughlen)];
        end
        if meanormedian==1
        roughmean(ii)=nanmean(roughtmp);
        else
        roughmean(ii)=nanmedian(roughtmp);
        end
    end
else
    sfacsprep=sfacspreall;
    sfacsjsp=sfacsjsall;
    roughmean=roughness_all;
end

rh_invjs=bspline_deboor(p+1,knots,sfacsjsp,t_rh);
rh_invpre=bspline_deboor(p+1,knots,sfacsprep,t_rh);

if numel(tgstring)>1
load(tgstring)
slvl=interp1(xaxis,slvl,t_rh,'linear');
tidey=slvl;
tidex=t_rh;
rh_invjs=nanmean(rh_invjs)-rh_invjs;
rh_invpre=nanmean(rh_invpre)-rh_invpre;
h_init=nanmean(h_init)-h_init;
in=isnan(tidey)==0 & isnan(rh_invjs)==0;
tidey=tidey-mean(tidey(in));
rms_js=rms(rh_invjs(in)-tidey(in));
in=isnan(tidey)==0 & isnan(rh_invpre)==0;
rms_pre=rms(rh_invpre(in)-tidey(in));
disp(['rms pre is ',num2str(rms_pre*100),' cm'])
disp(['rms js is ',num2str(rms_js*100),' cm'])
end

pointsperday=sum(~isnan(h_init))/(enddate+1-startdate);
disp(['points per day = ',num2str(pointsperday)])

if makefig==1
% PLOTTING
% Defaults
width =10;     % Width in inches % 3.5 was for putting two on one line i think
height =5;    % Height in inches
fsz = 13;      % Fontsize
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

figure('visible','on')
if roughnessplot==1
subplot(2,1,1)
end
if numel(tgstring)>0
%tideyt=tideyt-nanmean(tideyt);
scatter(tidex,tidey,1,'k')
hold on
end
scatter(x_init,h_init,'b+')
hold on
plot(t_rh,rh_invjs,'b','linewidth',2)
plot(t_rh,rh_invpre,'b--','linewidth',1)
axis([startdate enddate+1 -inf inf])
ylabel('Sea level (m)','interpreter','latex','fontsize',fsz)
datetick('x',1,'keeplimits','keepticks')
set(gca,'ticklabelinterpreter','latex','fontsize',fsz)
%h=legend('Tide gauge','GNSS-R: spectral analysis (unadjusted)','GNSS-R: inverse modelling');
%set(h,'interpreter','latex','fontsize',fsz)
if roughnessplot==1
subplot(2,1,2)
plot(startdate+0.5:1:enddate-0.5,roughmean,'r')
%plot(t_rh,movmean(rh_invjs(in)-tiderms(in),10),'r')
axis([startdate enddate+1 -inf inf])
datetick('x',1,'keeplimits','keepticks')
ylabel('Surface roughness_all (m)','interpreter','latex','fontsize',fsz)
set(gca,'ticklabelinterpreter','latex','fontsize',fsz)
%set(gca,'xtick',startdate:2:enddate)
end
print('invfig', '-dpng', '-r300');
end


end
