function [delay_out] = trop_delay_tp(curjd,lat,hell,hsfc,elv,p,Tm,e,ah,aw,lambda,psfc,Tmsfc,esfc)

% This function works out tropospheric delay using other functions that can
% be retrieved from https://vmf.geo.tuwien.ac.at/codes/

dmjd = curjd; % julian day
dlat(1) = lat*pi/180.d0;

zd = (90 - elv)*pi/180; % zenith distance in radians

% The zenith hydrostatic delay is always calculated with the equation by
% Saastamoinen (1972) as refined by Davis et al. (1985)
[zhd] = saasthyd (p,dlat,hell);
[zhdsfc] = saasthyd (psfc,dlat,hsfc);
delzhd=zhdsfc-zhd;

% The zenith wet delay can be calculated with (2) the equation 18 by Askne and
% Nordius (1987)
[zwd] = asknewet (e,Tm,lambda);
[zwdsfc] = asknewet (esfc,Tmsfc,lambda);
delzwd=zwdsfc-zwd;

% Calculate the mapping functions
[mfh,mfw] = vmf1_ht (ah,aw,dmjd,dlat,hell,zd);
        
% from Williams & Nievinski (2017)
delay_out=2*(delzhd*mfh + delzwd*mfw);

end
