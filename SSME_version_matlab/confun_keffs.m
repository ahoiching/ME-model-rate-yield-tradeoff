function [c, ceq] = confun_keffs(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis)
%CONFUN Summary of this function goes here
%   Detailed explanation goes here
mu = v(20)/sig;
%kres1=kres*1.3;
%kres1=kres;
kfer1=kfer;
%kfer1=kfer*1.3^(-mu);
% if mu<0.6
%     kfer1=kfer;
% else
%     kfer1=kfer/(mu*2);
% end

% if mu<0.77
%      kres1=kres;
% else
%      kres1=kres * mu * 1.3;
% end

kres1=kres*(1+1.3*mu^20)./(1+mu^20);
%kres1=kres;
kribo1=kribo;
%kfer1=kfer;
c = [
        v(2)*mu - v(14)*kres1
        v(3)*mu - v(15)*kfer1
        sum(v(7:12))*mu - v(16)*kribo1
    ];
%Nonlinear equality constraints
ceq = [
           phi0*mu*alpha*beta*(phi_max-lacZ) - v(18)
           mu*alpha*beta*(1-phi_max)-v(19)
           lacZ*mu*alpha*beta - v(17)
           sum(v(13:19)) - mu*beta
           sum(v(14:19)) - mu*alpha*beta
           eng_dis-v(21)
      ];
end

