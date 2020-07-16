function [A,G,names,ind] = coef2ampphase(coefs,names)

for ii=1:numel(coefs)/2
    A(ii)=sqrt(coefs(ii*2-1)^2+coefs(ii*2)^2);
    phi=mod(atan2(coefs(ii*2),coefs(ii*2-1))/pi*180,360); % y/x is y,x
    G(ii)=phi;
    report(ii,1)=cellstr([names(ii,:),' A=',num2str(round(A(ii)*100,4,'significant')),' g=',num2str(round(G(ii),4,'significant'))]);
end
[~,ind]=sort(A,'descend');
%A=A(ind);
%G=G(ind);
report=report(ind,:);
if numel(coefs)/2>=30
disp(report(1:30,1))
else
disp(report)
end
%names=names(ind,:);

end

