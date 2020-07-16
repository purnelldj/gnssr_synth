function [powerout] = scalenoisepower(station,powerin,freqsin)

% scale the noise (random sine waves) to match mean lsp shape
% need to add for each site

if strcmp(station,'GTGU')==1 %|| strcmp(station(1:4),'spby')
    p1 = -9.624e-10;
    p2 = 6.659e-07;
    p3 = -0.0001452;
    p4 = 0.008267;
    p5 = 0.5084;
    powerout=powerin.*polyval([p1 p2 p3 p4 p5],freqsin);
elseif strcmp(station(1:4),'sc02')==1
    powerout=powerin.*freqsin.^-0.8; % EXPONENTIAL
elseif strcmp(station(1:4),'bur2')==1
    powerout=powerin.*freqsin.^-0.3; % EXPONENTIAL
elseif strcmp(station(1:4),'scoa')==1
    powerout=powerin.*freqsin.^-0.3; % EXPONENTIAL
elseif strcmp(station(1:4),'spby')==1
    powerout=powerin.*freqsin.^-0.3; % EXPONENTIAL
elseif strcmp(station(1:4),'sab2')==1
    powerout=powerin.*freqsin.^-0.6; % EXPONENTIAL
else
    disp('you need to add the station to scalenoisepower.m ')
end

end

