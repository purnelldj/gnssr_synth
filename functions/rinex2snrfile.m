function [snr_data]=rinex2snrfile(obsstr,sp3str,elv_lims,azi_lims,staxyz,plat,plon)

%%

% THIS CODE TAKES RINEX FILES '.yyo' AND SATELLITE ORBIT DATA '.sp3'
% AND THEN OUTPUTS SOME ORAGNISED SNR DATA TO ANALYZE FOR REFLECTOMETRY
% WITH ELEVATION AND AZIMUTH ANGLE INPUTS

% INTERPOLATING SP3 FILES PROPERLY IS SO IMPORTANT
% I HAVE HAD SOME ISSUES WITH DIFFERENT SP3 FILES SO BE CAREFUL
% TRY USING DIFFERENT ONES FOR THE SAME TIME PERIOD AND THEN LOOK AT
% RESIDUALS

% NOTE
% THIS CODE WILL BE MUCH QUICKER IF YOU HAVE JUST S1 AND S2 DATA IN RINEX
% FILE, YOU CAN ACHIEVE THIS USING TEQC
% teqc +all -o.Obs s1+s2 FILE > FILE.YYo

s1data=NaN(86401,92);
s2data=NaN(86401,92);
satmatr={'G01';'G02';'G03';'G04';'G05';'G06';'G07';'G08';'G09';...
    'G10';'G11';'G12';'G13';'G14';'G15';'G16';'G17';'G18';'G19';...
    'G20';'G21';'G22';'G23';'G24';'G25';'G26';'G27';'G28';'G29';...
    'G30';'G31';'G32';...
    'R01';'R02';'R03';'R04';'R05';'R06';'R07';'R08';'R09';...
    'R10';'R11';'R12';'R13';'R14';'R15';'R16';'R17';'R18';'R19';...
    'R20';'R21';'R22';'R23';'R24';...
    'E01';'E02';'E03';'E04';'E05';'E06';'E07';'E08';'E09';...
    'E10';'E11';'E12';'E13';'E14';'E15';'E16';'E17';'E18';'E19';...
    'E20';'E21';'E22';'E23';'E24';'E25';'E26';'E27';'E28';'E29';...
    'E30';'E31';'E32';'E33';'E34';'E35';'E36'}; 

fid=fopen(obsstr);

tline=fgets(fid);
endnow=0;
disp('reading rinex header')
while endnow==0 % find end of header
    tline=fgets(fid);
    obsq=strfind(tline,'# / TYPES OF OBSERV');
    if size(obsq,1)>0
        numobs=str2double(tline(1:6));
        if numobs<6
            obsl=1;
        elseif numobs>5 && numobs<11
            obsl=2;
        elseif numobs>10 && numobs<16
            obsl=3;
        elseif numobs>15 && numobs<21
            obsl=4;
        elseif numobs>20 && numobs<26
            obsl=5;
        end
        clear s1pos s2pos
        for ii=1:numobs
            if ii==10 || ii==19
                tline=fgets(fid);
            end
            
            if ii<10
            obstype(ii,:)=tline(ii*6+5:ii*6+6);
            elseif ii>9 && ii<19
            obstype(ii,:)=tline((ii-9)*6+5:(ii-9)*6+6);
            elseif ii>18
            obstype(ii,:)=tline((ii-18)*6+5:(ii-18)*6+6);
            end
            if strcmp(obstype(ii,:),'S1')==1
                s1pos=ii;
            elseif strcmp(obstype(ii,:),'S2')==1
                s2pos=ii;
            end
        end
    end
    endq=strfind(tline,'END OF HEADER');
    if size(endq,1)>0
        endnow=1;
        tline=fgets(fid);
    end
end

% NOW SORTING OUT THE OBSERVATIONS
disp('extracting snr data')
cursecs=0;
numsats=0;
issue=0;
while ~feof(fid)
    i=0;
    % this is to deal with spliced files
    if numel(tline)>=31
    if strcmp(tline(1:31),'                            4  ')==1 ||...
            strcmp(tline(1:32),'                            4 64')
        tmpskip=str2double(tline(32))+1;
        while i<tmpskip
            tline=fgets(fid);
            i=i+1;
        end
        if numel(tline)>=67
        while strcmp(tline(61:67),'COMMENT')
            tline=fgets(fid);
            i=i+1;
        end
        end
        continue
    end
    end

    cursecsold=cursecs;
    cursecs=str2double(tline(11:12))*60*60+str2double(tline(14:15))*60+str2double(tline(16:26));
    numsatsold=numsats;
    numsats=str2double(tline(31:32));
    if isnan(numsats)==1 || isnan(cursecs)==1 || (cursecs==0 && cursecsold~=0)
        issue=1;
        disp('issue')
        return
    end
    if numel(tline)<34
        while i<numsats+1
            tline=fgets(fid);
            i=i+1;
        end
        continue
    end
    if cursecs>86400 %|| cursecs-cursecsold>1
        break
    end
    % here just saving the order of the satellites
    numsatst=min([12 numsats]);
    satsord=tline(33:32+numsatst*3);
    mult12=12;
    ab=1;
    while mult12<numsats
        ab=ab+1;
        mult12=12*ab;
        numsatst=min([mult12 numsats])-12*(ab-1);
        tline=fgets(fid);
        satsord=[satsord,tline(33:32+(numsatst)*3)];
    end
    tline=fgets(fid);
    
    % now collect the data
    for ss=1:numsats
        satname=satsord((ss-1)*3+1:(ss-1)*3+3);
        satis=strcmp(satname,satmatr(:,1))==1;
        satind=sum(satis.*[1:size(satis,1)].');
        i=0;
        if sum(satis)==0
           while i~=obsl
             tline=fgets(fid);
             i=i+1;
            end
            continue
        end
        % getting s1
        mult5=5;
        while s1pos>mult5
            i=i+1;
            mult5=5*(i+1);
            tline=fgets(fid);
        end
        if numel(tline)>=(s1pos-i*5-1)*16+14
        s1data(cursecs+1,satind)=str2double(tline((s1pos-i*5-1)*16+1:(s1pos-i*5-1)*16+14));
        end
        % getting s2
        while s2pos>mult5
            i=i+1;
            mult5=5*(i+1);
            tline=fgets(fid);
        end
        if numel(tline)>=(s2pos-i*5-1)*16+14
        s2data(cursecs+1,satind)=str2double(tline((s2pos-i*5-1)*16+1:(s2pos-i*5-1)*16+14));
        end
        while i~=obsl
            tline=fgets(fid);
            i=i+1;
        end
    end
end

disp('getting orbit info')

[txyz,xyz] = readsp3file(sp3str);
txyzsecs=txyz-txyz(1);
txyzsecs=txyzsecs.*86400;

snr_data=[];

for satind=1:92

s1datat=squeeze(s1data(:,satind));
s2datat=squeeze(s2data(:,satind));
secs=find(~isnan(s1datat(:)) | ~isnan(s2datat(:)));
s1datat=s1datat(secs);
s2datat=s2datat(secs);

clear xyzt
if sum(~isnan(squeeze(xyz(satind,:,1))))>1
xyzt(:,1)=spline(txyzsecs,squeeze(xyz(satind,:,1)),secs-1);
xyzt(:,2)=spline(txyzsecs,squeeze(xyz(satind,:,2)),secs-1);
xyzt(:,3)=spline(txyzsecs,squeeze(xyz(satind,:,3)),secs-1);
else
    continue
end

ind=0;
snrdatat=[];
for tt=1:numel(secs)
[azi,elv]=gnss2azelv(staxyz,xyzt(tt,:),plat,plon);
if (azi>azi_lims(1) && azi<azi_lims(2)) && (elv>elv_lims(1) && elv<elv_lims(2))
ind=ind+1;
snrdatat(ind,1)=satind;
snrdatat(ind,2)=elv;
snrdatat(ind,3)=azi;
snrdatat(ind,4)=secs(tt)-1;
snrdatat(ind,5)=NaN;
snrdatat(ind,6)=NaN;
snrdatat(ind,7)=s1datat(tt);
snrdatat(ind,8)=s2datat(tt);
end
end

snr_data=[snr_data;snrdatat];

end

if size(snr_data,1)>0
snr_data=sortrows(snr_data,4);
else
    disp('sorry no data mate')
end

disp('*************done*************')

end


