function tidesout = tidemod_kl_plot(coefs,t,rh_nadj,ju,tanthter,freqs,lat)

% times should be in datenum format (days)
times0=t-median(t);

ltype='nodal';
ctime=median(t);
%addpath('t_tide_v1.4beta/')
[v,u,f]=t_vuf(ltype,ctime,ju+1,lat);
v=v.*360;
u=u.*360;

tidesout=zeros(numel(t),1);
for ii=1:numel(freqs)
    tidesout=tidesout+2*pi*freqs(ii)*f(ii).*(-coefs(ii*2-1)*sind(360*freqs(ii).*times0.*24+u(ii)+v(ii))...
        +coefs(ii*2)*cosd(360*freqs(ii).*times0.*24+u(ii)+v(ii))).*tanthter;
end
tidesout=rh_nadj-tidesout;

end
