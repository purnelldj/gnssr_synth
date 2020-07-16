% this code extracts GNSS satellite orbit data from nasa archive

clear all

% change this to wherever your 'gnss functions' directory is
%addpath('/codes/gnss_functions')

station='scoa';
datevecs=datenum(2019,01,01);
datevece=datenum(2019,01,02);
curdatevec=datevecs-1;

while curdatevec<datevece
    curdatevec=curdatevec+1;
    curdt=datetime(curdatevec,'convertfrom','datenum');
    strday=char(datetime(curdt,'format','DDD'));
    stryrs=char(datetime(curdt,'format','yy'));
    stryrl=char(datetime(curdt,'format','yyyy'));
    ftpobj = ftp('rgpdata.ign.fr');
    cd(ftpobj,['pub/data/',stryrl,'/',strday,'/data_1']);
    mget(ftpobj,[station,strday,'*.',stryrs,'d.Z']);
    close(ftpobj)
end

%%




