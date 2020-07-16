% this code extracts GNSS satellite orbit data from nasa archive

clear all

% change this to wherever your 'gnss functions' directory is
%addpath('/codes/gnss_functions')

station='bur2';
datevecs=datenum(2019,05,01);
datevece=datenum(2019,06,01);
curdatevec=datevecs-1;
hhdir=['00';'01';'02';'03';'04';'05';'06';'07';'08';'09';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22';'23'];

while curdatevec<datevece
    curdatevec=curdatevec+1;
    curdt=datetime(curdatevec,'convertfrom','datenum');
    disp(char(curdt))
    strday=char(datetime(curdt,'format','DDD'));
    stryrs=char(datetime(curdt,'format','yy'));
    stryrl=char(datetime(curdt,'format','yyyy'));
    for hh=1:24
    ftpobj = ftp('ftp.ga.gov.au');
    cd(ftpobj,['geodesy-outgoing/gnss/data/highrate/',stryrl,'/',stryrs,strday,'/',hhdir(hh,:)]);
    mget(ftpobj,[station,strday,'*.',stryrs,'d.Z']);
    close(ftpobj)
    end
end

return

%% to get the nav data

station='bur2';
datevecs=datenum(2019,01,1);
datevece=datenum(2019,02,01);
curdatevec=datevecs-1;

while curdatevec<datevece
    curdatevec=curdatevec+1;
    curdt=datetime(curdatevec,'convertfrom','datenum');
    disp(char(curdt))
    strday=char(datetime(curdt,'format','DDD'));
    stryrs=char(datetime(curdt,'format','yy'));
    stryrl=char(datetime(curdt,'format','yyyy'));
    ftpobj = ftp('ftp.ga.gov.au');
    cd(ftpobj,['geodesy-outgoing/gnss/data/daily/',stryrl,'/',stryrs,strday,'/']);
    mget(ftpobj,[station,strday,'*.',stryrs,'n.Z']);
    close(ftpobj)
end




