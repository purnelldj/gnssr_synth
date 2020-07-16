function [res] = bspline_js(coefs_0,t1_all,sinelv1_all,snr1_all,knots,bspline_order,...
    satno_all,gps,glo,gal,antno_all,meanhgts)

t1_alls=t1_all;
sinelv1_alls=sinelv1_all;
snr1_alls=snr1_all;
satno_alls=satno_all;
%antno_alls=antno_all;

res=[];
ants=unique(antno_all);

consts=gps+glo+gal;
height_coeff=coefs_0(1:end-consts*2-1);

for k=1:numel(ants)

now=antno_all(:)==k;
t1_all=t1_alls(now);
sinelv1_all=sinelv1_alls(now);
snr1_all=snr1_alls(now);
satno_all=satno_alls(now);

h1=bspline_deboor(bspline_order+1,knots,height_coeff,t1_all);
h1=h1.';
if k>1
    h1=h1+meanhgts(k)-meanhgts(1);
end

tmpc=0;
% GPS
if gps==1
tmpc=tmpc+1;
gps=satno_all(:,1)<33;
L1car=(299792458/(1575.42e06/1.023)); % for GPS
L1k=(2*pi)/L1car;
modelSNR1 = (coefs_0(end-1)*sin((4*pi*h1(gps).*sinelv1_all(gps))/L1car)+...
    coefs_0(end-2)*cos((4*pi*h1(gps).*sinelv1_all(gps))/L1car)).*...
    exp(-4*L1k^2*coefs_0(end)*sinelv1_all(gps).^2);
res=[res;modelSNR1-snr1_all(gps)];
end
% GLO
if glo==1
tmpc=tmpc+1;
glo=satno_all(:,1)>32 & satno_all(:,1)<57;
load('glonasswlen.mat')
satnotmp=satno_all(glo,1);
for ij=1:numel(satnotmp)
L1car(ij)=glonasswlen(satnotmp(ij)-32);
end
L1car=L1car.';
L1k=(2*pi)./L1car;
modelSNR1 = (coefs_0(end-(tmpc-1)*2-1)*sin((4*pi*h1(glo).*sinelv1_all(glo))./L1car)+...
    coefs_0(end-(tmpc-1)*2-2)*cos((4*pi*h1(glo).*sinelv1_all(glo))./L1car)).*...
    exp(-4*L1k.^2*coefs_0(end).*sinelv1_all(glo).^2);
res=[res;modelSNR1-snr1_all(glo)];
end
% GAL
if gal==1
tmpc=tmpc+1;
gal=satno_all(:,1)>56;
L1car=(299792458/(1575.42e06/1.023)); % for Galileo
L1k=(2*pi)/L1car;
modelSNR1 = (coefs_0(end-(tmpc-1)*2-1)*sin((4*pi*h1(gal).*sinelv1_all(gal))/L1car)+...
    coefs_0(end-(tmpc-1)*2-2)*cos((4*pi*h1(gal).*sinelv1_all(gal))/L1car)).*...
    exp(-4*L1k^2*coefs_0(end)*sinelv1_all(gal).^2);
res=[res;modelSNR1-snr1_all(gal)];
end
end

end

