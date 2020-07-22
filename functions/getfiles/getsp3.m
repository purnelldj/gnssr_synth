% this code extracts GNSS satellite orbit data from nasa archive

clear all

% inputs
gpsweeks=2108; % will download from this week
gpsweeke=2114; % to this week

gpsweekn=gpsweeks-1;
while gpsweekn<gpsweeke
    gpsweekn=gpsweekn+1;
    gpsweek=num2str(gpsweekn);

    disp(['gpsweekn=',gpsweek])
    
% use this script for downloading sp3 data
if gpsweekn < 1962
ftpobj = ftp('cddis.gsfc.nasa.gov');
cd(ftpobj,strcat('pub/gps/products/mgex/',gpsweek));
dow=0;
while dow<7
mget(ftpobj,strcat('com',gpsweek,num2str(dow),'.sp3.Z'));
dow=dow+1;
end
else
%ftpobj = ftp('cddis.gsfc.nasa.gov/pub/gps/products/mgex/');
ftpobj = ftp('cddis.gsfc.nasa.gov');
cd(ftpobj,strcat('pub/gps/products/mgex/',gpsweek));
%cd(ftpobj,strcat('pub/gnss/products/',gpsweek));
dow=0;
sow=0;
while dow<7
jd=gps2jd(gpsweekn,sow,0);
curdt=datetime(jd,'convertfrom','juliandate');
strday=char(datetime(curdt,'format','DDD'));
stryrl=char(datetime(curdt,'format','yyyy'));
mget(ftpobj,strcat('COD0MGXFIN_',stryrl,strday,'0000_01D_05M_ORB.SP3.gz'));
%mget(ftpobj,strcat('GFZ0MGXRAP_',year,doys,'0000_01D_05M_ORB.SP3.gz'));
%mget(ftpobj,strcat('igs',gpsweek,num2str(dow),'.sp3.Z'));
dow=dow+1;
sow=dow*86400;
end
end
end

