function resid = noisefunresid(coefs,snr_preproc,powerin,freqsin,sinelv,phaseoff,resln,pkin,varint,maxf1,ovs)
snr_recr=noisefun(coefs,snr_preproc,powerin,freqsin,sinelv,phaseoff,resln);
fi=1:1:maxf1;
[psdn,~,~,~,~,~]=fLSPw(sinelv,snr_recr,fi,0.05,ovs);
psdn=cell2mat(psdn);
varresid=var(snr_recr)-varint;
pkresid=max(psdn)-pkin;
resid=[varresid;pkresid];
end
