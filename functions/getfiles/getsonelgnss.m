% this code extracts GNSS satellite orbit data from nasa archive

clear all

% change this to wherever your 'gnss functions' directory is
%addpath('/codes/gnss_functions')

station='scoa';
datevecs=datenum(2019,05,01);
datevece=datenum(2019,06,01);
curdatevec=datevecs-1;

while curdatevec<datevece
    curdatevec=curdatevec+1;
    curdt=datetime(curdatevec,'convertfrom','datenum');
    strday=char(datetime(curdt,'format','DDD'));
    stryrs=char(datetime(curdt,'format','yy'));
    stryrl=char(datetime(curdt,'format','yyyy'));
    ftpobj = ftp('ftp.sonel.org');
    cd(ftpobj,['gps/data/',stryrl,'/',strday]);
    mget(ftpobj,[station,strday,'0.',stryrs,'d.Z']);
    close(ftpobj)
end

