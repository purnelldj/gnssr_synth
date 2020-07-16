% this code extracts GNSS satellite orbit data from nasa archive

clear all

% change this to wherever your 'gnss functions' directory is
addpath('/Users/dave/Library/Mobile Documents/com~apple~CloudDocs/work/gnssr_matlab_NEW/functions/geodetic299/geodetic')

% inputs
gpsweeks=2096; % will download from this week
gpsweeke=2108; % to this week

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
    year=jd2yr(jd);
    year=num2str(year);
    year=year(1:4);
    doy=jd2doy(jd);
    if doy<100
        if doy<10
            doys=strcat('00',num2str(doy));
        else
            doys=strcat('0',num2str(doy));
        end
    else
        doys=num2str(doy);
    end
mget(ftpobj,strcat('COD0MGXFIN_',year,doys,'0000_01D_05M_ORB.SP3.gz'));
%mget(ftpobj,strcat('GFZ0MGXRAP_',year,doys,'0000_01D_05M_ORB.SP3.gz'));
%mget(ftpobj,strcat('igs',gpsweek,num2str(dow),'.sp3.Z'));
dow=dow+1;
sow=dow*86400;
end
end
end

