function [snr_recr2,yt2] = noisefun(coefs,snr_preproc,powerin,freqsin,sinelv,phaseoff,resln)

coefs(1)=exp(coefs(1));
coefs(2)=exp(coefs(2));

SNRt=sqrt(10.^(snr_preproc./10));
pt=polyfit(sinelv,SNRt,2);
yt=polyval(pt,sinelv);
SNRdt=SNRt-yt;
powerint=coefs(1).*powerin;
snr_recr=SNRdt*coefs(2);
for ii=1:numel(freqsin)
    snr_recr=snr_recr+powerint(ii)*sin(freqsin(ii)*(2*pi).*sinelv+phaseoff(ii));
end
snr_recr=snr_recr-mean(snr_recr);
snr_recr2=snr_recr+yt;
snr_recr2=10.*log10(snr_recr2.^2);
reslnon=0;
if reslnon==1
mods=mod(snr_recr2,resln);
for i=1:size(mods,2)
    if mods(i)>=resln/2
        snr_recr2(i)=snr_recr2(i)-mods(i)+resln;
    else
        snr_recr2(i)=snr_recr2(i)-mods(i);
    end
end
end
snr_recr2=sqrt(10.^(snr_recr2./10));
pt2=polyfit(sinelv,snr_recr2,2);
yt2=polyval(pt2,sinelv);
snr_recr2=snr_recr2-yt2;

end

