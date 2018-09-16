function [c, ceq] = confun(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis)
%CONFUN Summary of this function goes here
%   Detailed explanation goes here
mu = v(20)/sig;
c = [
        v(2)*mu - v(14)*kres
        v(3)*mu - v(15)*kfer
        sum(v(7:12))*mu - v(16)*kribo
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

