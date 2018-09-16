clc;
%%
%data from NQ381,and NCM3722
CUF_NQ381=[28.56,36.72,39.24]; %carbon uptake flux
MUL_NQ381=[0.77,0.82,0.91]; %growth rate
ACT_NQ381=[1.32/3,4.05/3,9.18/3]; %acetate secretion
CO2_r_NQ381=[8.49,6.90,4.90];%co2 secretion from respirotory

CUF_NCM3722=[41.52]; %carbon uptake flux
MUL_NCM3722=[0.97]; %growth rate
ACT_NCM3722=[9.72/3]; %acetate secretion
CO2_r_NCM3722=[3.47];%co2 secretion from respirotory
%Part 1: show that this small-scale ME-model creates equivalent result as Basan's model
%Parameters from Basan et al. 2015
beta=28.5;
alpha=1;
eps_r=390;
eps_f=750;
e_r=4.4;
e_f=2;
b=0.12;
sig=45.7;
phi0=0.81;
phi_max=0.42;
lacZ=0;
eng_dis=0;

x0 = double(ones(21,1)); %Make a staring guess at the solution
options = optimoptions (@fmincon, 'Algorithm', 'sqp');
S_struct = load('S(basan).mat');
Aeq=S_struct.S;
beq=double(zeros(1,11));
lb = double(zeros(21,1));
ub = double(1000 * ones(21,1));

loopcount=500;
mul=zeros(1,loopcount);
act_ex=zeros(1,loopcount);
co2_r_ex=zeros(1,loopcount);

surl = zeros(1,loopcount);
sur = 1;
increment = 0.2;

%The derivation of keffs in ME-model
kres=eps_r/(e_r*beta*alpha*(phi_max));
kfer=eps_f/(e_f*beta*alpha*(phi_max));
kribo=1/(b*(phi_max));


nonlcon = @(v)confun(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis);

for index=1:loopcount
   surl(index)=sur;
   ub(1)=sur;
   exitflag=0;
   while exitflag ~= 1
       x0 = x0*1.01; 
       [v,fval,exitflag]=fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub,nonlcon,options);
   end
   mul(index)= v(20)/sig;
   co2_r_ex(index)=v(4);
   act_ex(index)=v(5);
   sur=sur+increment;
   x0 = v;
end

%sampling the growth rate list and make it sparser
mul_sparse=mul(1:10:500);
%parameters in basan's model:
S_ac=1/3;
S_co2_r=1;
%data from basan's model:
mu_co2=(1-phi0)/(sig/eps_f+b);
mu_act=(1-phi0)/(sig/eps_r+b);
co2_r_ex_d=(mul_sparse<=mu_co2).*(S_co2_r*(eps_r/e_r)*((1-phi0)-(sig/eps_f+b)*mul_sparse)/(1-(eps_r/eps_f)));
act_ex_d=(mul_sparse>=mu_act).*(S_ac*(eps_f/e_f)*((sig/eps_r+b)*mul_sparse-(1-phi0))/((eps_f/eps_r)-1));
figure;
plot(mul,co2_r_ex,'b',mul,act_ex,'r',mul_sparse,co2_r_ex_d,'b.',mul_sparse,act_ex_d,'r.',...
     MUL_NQ381,ACT_NQ381,'ro',MUL_NQ381,CO2_r_NQ381,'bo',...
     MUL_NCM3722,ACT_NCM3722,'r^',MUL_NCM3722,CO2_r_NCM3722,'b^',...
     'MarkerSize',10);
title('Growth-rate dependence of acetate excretion and respiratory CO_2 production(chemostat)');
xlabel('Growth rate, \mu(h^{-1})');
ylabel('Flux (mmol OD600^{-1} h^{-1})');
legend('ME:v_{EX\_CO_2,r}','ME:v_{EX\_acetate}',...
       'Basan:v_{EX\_CO_2,r}','Basan:v_{EX\_acetate}',...
       'NQ381 data:v_{EX\_acetate}','NQ381 data:v_{EX\_CO_2,r}',...
       'NCM3722 data:v_{EX\_acetate}','NCM3722 data:v_{EX\_CO_2,r}',...
       'Location','northeast');
   
%%
%Part 2: The effect of alpha:
loopcount=50;
muary=zeros(loopcount,5);
co2ary=zeros(loopcount,5);
actary=zeros(loopcount,5);
proary=zeros(loopcount,5);
obmsary=zeros(loopcount,5);
for outloop=1:5
   alpha=outloop/5; 
   mul=zeros(1,loopcount);
   %engdil=zeros(1,loopcount);
   act_ex=zeros(1,loopcount);
   co2_r_ex=zeros(1,loopcount);
   prot_dmd=zeros(1,loopcount);
   obms_dmd=zeros(1,loopcount);
   %rp_ratio=zeros(1,loopcount);

   surl = zeros(1,loopcount);
   sur = 1;
   increment = 2;

   kres=eps_r/(e_r*beta*alpha*(phi_max));
   kfer=eps_f/(e_f*beta*alpha*(phi_max));
   kribo=1/(b*(phi_max));

   nonlcon = @(v)confun(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis);
   for index=1:loopcount
     surl(index)=sur;
     ub(1)=sur;
     exitflag=0;
     while exitflag ~= 1
         x0 = x0*1.01; 
         [v,fval,exitflag]=fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub,nonlcon,options);
     end
     mul(index)= v(20)/sig;
     co2_r_ex(index)=v(4);
     act_ex(index)=v(5);
     prot_dmd(index)=sum(v(14:19));
     obms_dmd(index)=v(13);
     sur=sur+increment;
     x0 = v;
   end
   muary(1:loopcount,outloop)=mul;
   co2ary(1:loopcount,outloop)=co2_r_ex;
   actary(1:loopcount,outloop)=act_ex;
   proary(1:loopcount,outloop)=prot_dmd;
   obmsary(1:loopcount,outloop)=obms_dmd;
end
figure;
plot(surl,muary,CUF_NQ381,MUL_NQ381,'S',CUF_NCM3722,MUL_NCM3722,'o','LineWidth',3);
title({'Growth rates on varied substrate uptake bound','with different proteome/biomass(\alpha) values'});
xlabel('Substrate(Carbon) uptake rate bound(mmol OD600^{-1} h^{-1})');
ylabel('Growth rate, \mu(h^{-1})');
legend('\alpha=0.2','\alpha=0.4','\alpha=0.6','\alpha=0.8','\alpha=1','NQ381 data','NCM3722 data','Location','northwest');

figure;
plot(muary,co2ary,'.',muary,actary,'o',muary,proary,'x--',muary,obmsary,'^');
title({'Growth-rate dependence of acetate excretion, respiratory CO_2 production, as well as total protein production','with different proteome/biomass(\alpha) values'});
xlabel('Growth rate, \mu(h^{-1})');
ylabel('Flux (mmol OD600^{-1} h^{-1})');
legend('CO_2,\alpha=0.2','CO_2,\alpha=0.4','CO_2,\alpha=0.6','CO_2,\alpha=0.8','CO_2,\alpha=1',...
       'Act,\alpha=0.2','Act,\alpha=0.4','Act,\alpha=0.6','Act,\alpha=0.8','Act,\alpha=1',...
       'pro,\alpha=0.2','pro,\alpha=0.4','pro,\alpha=0.6','pro,\alpha=0.8','pro,\alpha=1',...
       'other bms,\alpha=0.2','other bms,\alpha=0.4','other bms,\alpha=0.6','other bms,\alpha=0.8','other bms,\alpha=1',...
       'Location','northeast');

%%
%Part3: variable sensitivity, keffs:
alpha=0.5;
%(1) kres:
kres0=eps_r/(e_r*beta*alpha*(phi_max));
kfer=eps_f/(e_f*beta*alpha*(phi_max));
kribo=1/(b*(phi_max));
loopcount=500;

kresl=zeros(1,5);
for outloop=1:5
   kres=kres0/(1.2^(outloop-2)); 
   kresl(outloop)=kres;
   mul=zeros(1,loopcount);
   %engdil=zeros(1,loopcount);
   act_ex=zeros(1,loopcount);
   co2_r_ex=zeros(1,loopcount);
   %prot_dmd=zeros(1,loopcount);
   %rp_ratio=zeros(1,loopcount);

   surl = zeros(1,loopcount);
   sur = 1;
   increment = 0.2;

   nonlcon = @(v)confun(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis);
   for index=1:loopcount
   surl(index)=sur;
   ub(1)=sur;
   [v,fval]=fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub,nonlcon,options);
   mul(index)= v(20)/sig;
   co2_r_ex(index)=v(4);
   act_ex(index)=v(5);
   %prot_dmd(index)=sum(v(11:15));
   sur=sur+increment;
   x0 = v;
   end
   muary(1:loopcount,outloop)=mul;
   co2ary(1:loopcount,outloop)=co2_r_ex;
   actary(1:loopcount,outloop)=act_ex;
   %proary(1:loopcount,outloop)=prot_dmd;
end

figure;
h=zeros(1,5);
colors =colormap(hsv(5));
hold on
for k=1:5
    h(k)=plot(surl,muary(1:loopcount,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('k_{eff,res}=%f',kresl(k)));
end
title({'Growth rates on varied substrate uptake bound','with different k_{eff,res} values'});
xlabel('Substrate(Carbon) uptake rate bound(mmol OD600^{-1} h^{-1})');
ylabel('Growth rate, \mu(h^{-1})');
legend(h,'Location','NorthWest');
hold off;
figure;
%h=zeros(3,5);
h=zeros(2,5);
colors =colormap(hsv(5));
hold on
for k=1:5
    h(1,k)=plot(muary(1:loopcount,k),co2ary(1:loopcount,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('co2:k_{eff,res}=%f',kresl(k)));
    h(2,k)=plot(muary(1:loopcount,k),actary(1:loopcount,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('act:k_{eff,res}=%f',kresl(k)));
%     h(3,k)=plot(muary(1:loopcount,k),proary(1:loopcount,k),'o-','Color',colors(k,:),...
%         'DisplayName',sprintf('pro:kres=%f',kresl(k)));
end
legend(h(:),'Location','NorthWest');
%title({'Growth-rate dependence of acetate excretion, respiratory CO_2 production, as well as total protein production','with different kres values'});
title({'Growth-rate dependence of acetate excretion and respiratory CO_2 production','with different k_{eff,res} values'});
xlabel('Growth rate, \mu(h^{-1})');
ylabel('Flux (mmol OD600^{-1} h^{-1})');
hold off;

%(2) kfer:
kres=eps_r/(e_r*beta*alpha*(phi_max));
kfer0=eps_f/(e_f*beta*alpha*(phi_max));
kribo=1/(b*(phi_max));
loopcount=500;

kferl=zeros(1,5);
for outloop=1:5
   kfer=kfer0/(1.2^(outloop-2)); 
   kferl(outloop)=kfer;
   mul=zeros(1,loopcount);
   %engdil=zeros(1,loopcount);
   act_ex=zeros(1,loopcount);
   co2_r_ex=zeros(1,loopcount);
   %prot_dmd=zeros(1,loopcount);
   %rp_ratio=zeros(1,loopcount);

   surl = zeros(1,loopcount);
   sur = 1;
   increment = 0.2;

   nonlcon = @(v)confun(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis);
   for index=1:loopcount
   surl(index)=sur;
   ub(1)=sur;
   [v,fval]=fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub,nonlcon,options);
   mul(index)= v(20)/sig;
   co2_r_ex(index)=v(4);
   act_ex(index)=v(5);
   %prot_dmd(index)=sum(v(11:15));
   sur=sur+increment;
   x0 = v;
   end
   muary(1:loopcount,outloop)=mul;
   co2ary(1:loopcount,outloop)=co2_r_ex;
   actary(1:loopcount,outloop)=act_ex;
   %proary(1:loopcount,outloop)=prot_dmd;
end

figure;
h=zeros(1,5);
colors =colormap(hsv(5));
hold on
for k=1:5
    h(k)=plot(surl,muary(1:loopcount,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('k_{eff,fer}=%f',kferl(k)));
end
title({'Growth rates on varied substrate uptake bound','with different k_{eff,fer} values'});
xlabel('Substrate(Carbon) uptake rate bound(mmol OD600^{-1} h^{-1})');
ylabel('Growth rate, \mu(h^{-1})');
legend(h,'Location','NorthWest');
hold off;
figure;
%h=zeros(3,5);
h=zeros(2,5);
colors =colormap(hsv(5));
hold on
for k=1:5
    h(1,k)=plot(muary(1:loopcount,k),co2ary(1:loopcount,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('co2:k_{eff,fer}=%f',kferl(k)));
    h(2,k)=plot(muary(1:loopcount,k),actary(1:loopcount,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('act:k_{eff,fer}=%f',kferl(k)));
%     h(3,k)=plot(muary(1:loopcount,k),proary(1:loopcount,k),'o-','Color',colors(k,:),...
%         'DisplayName',sprintf('pro:kfer=%f',kferl(k)));
end
legend(h(:),'Location','NorthWest');
title({'Growth-rate dependence of acetate excretion and respiratory CO_2 production','with different k_{eff,fer} values'});
xlabel('Growth rate, \mu(h^{-1})');
ylabel('Flux (mmol OD600^{-1} h^{-1})');
hold off;

%(3) kribo(translation limitation);
kres=eps_r/(e_r*beta*alpha*(phi_max));
kfer=eps_f/(e_f*beta*alpha*(phi_max));
kribo0=1/(b*(phi_max));
loopcount=500;

kribol=zeros(1,5);
for outloop=1:5
   kribo=kribo0/(1.2^(outloop-2)); 
   kribol(outloop)=kribo;
   mul=zeros(1,loopcount);
   %engdil=zeros(1,loopcount);
   act_ex=zeros(1,loopcount);
   co2_r_ex=zeros(1,loopcount);
   %prot_dmd=zeros(1,loopcount);
   %rp_ratio=zeros(1,loopcount);

   surl = zeros(1,loopcount);
   sur = 1;
   increment = 0.2;

   nonlcon = @(v)confun(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis);
   for index=1:loopcount
   surl(index)=sur;
   ub(1)=sur;
   [v,fval]=fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub,nonlcon,options);
   mul(index)= v(20)/sig;
   co2_r_ex(index)=v(4);
   act_ex(index)=v(5);
   %prot_dmd(index)=sum(v(11:15));
   sur=sur+increment;
   x0 = v;
   end
   muary(1:loopcount,outloop)=mul;
   co2ary(1:loopcount,outloop)=co2_r_ex;
   actary(1:loopcount,outloop)=act_ex;
   %proary(1:loopcount,outloop)=prot_dmd;
end

figure;
h=zeros(1,5);
colors =colormap(hsv(5));
hold on
for k=1:5
    h(k)=plot(surl,muary(1:loopcount,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('k_{eff,ribo}=%f',kribol(k)));
end
title({'Growth rates on varied substrate uptake bound','with different k_{eff,ribo} values'});
xlabel('Substrate(Carbon) uptake rate bound(mmol OD600^{-1} h^{-1})');
ylabel('Growth rate, \mu(h^{-1})');
legend(h,'Location','NorthWest');
hold off;
figure;
% h=zeros(3,5);
h=zeros(2,5);
colors =colormap(hsv(5));
hold on
for k=1:5
    h(1,k)=plot(muary(1:loopcount,k),co2ary(1:loopcount,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('co2:k_{eff,ribo}=%f',kribol(k)));
    h(2,k)=plot(muary(1:loopcount,k),actary(1:loopcount,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('act:k_{eff,ribo}=%f',kribol(k)));
%     h(3,k)=plot(muary(1:loopcount,k),proary(1:loopcount,k),'o-','Color',colors(k,:),...
%         'DisplayName',sprintf('pro:kribo=%f',kribol(k)));
end
legend(h(:),'Location','NorthWest');
title({'Growth-rate dependence of acetate excretion, respiratory CO_2 production, as well as total protein production','with different k_{eff,ribo} values'});
xlabel('Growth rate, \mu(h^{-1})');
ylabel('Flux (mmol OD600^{-1} h^{-1})');
hold off;

%%
%Part 4: Recreation of the graph in the paper.
% beta=28.5;
% alpha=1;
% eps_r=260;
% eps_f=2600;
% e_r=4.4;
% e_f=2;
% b=0.055;
% sig=45.7;
% phi0=0.81;
% phi_max=0.42;
% lacZ=0;
% eng_dis=0;

%(1)LacZ
MUL_NQ1389_phiZ_8=[0.93 0.8 0.7 0.73];
ACT_NQ1389_phiZ_8=[3.01 1.33 0.68 0.7];
MUL_NQ1389_phiZ_16=[0.72 0.61 0.55 0.51];
ACT_NQ1389_phiZ_16=[2.35 1.02 0.53 0.43];

beta=28.5;
alpha=1;
% eps_r=390;
% eps_f=750;
% e_r=4.4;
% e_f=2;
% b=0.13;
sig=45.7;
phi0=0.81;
phi_max=0.42;
lacz0=0;
eng_dis=0;
loopcount=500;

laczl=zeros(1,5);
for outloop=1:5
   lacZ=lacz0+(outloop-1)*0.04; 
   %lacZ=lacz0; 
   laczl(outloop)=lacZ;
   mul=zeros(1,loopcount);
   %engdil=zeros(1,loopcount);
   act_ex=zeros(1,loopcount);
   co2_r_ex=zeros(1,loopcount);
   %prot_dmd=zeros(1,loopcount);
   %rp_ratio=zeros(1,loopcount);

   surl = zeros(1,loopcount);
   sur = 1;
   increment = 0.2;
   
%    kres=eps_r/(e_r*beta*alpha*(phi_max));
%    kfer=eps_f/(e_f*beta*alpha*(phi_max));
%    kribo=1/(b*(phi_max));
   
   kres=3.8;
   kfer=41;
   kribo=195;

   nonlcon = @(v)confun(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis);
   for index=1:loopcount
   surl(index)=sur;
   ub(1)=sur;
   [v,fval]=fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub,nonlcon,options);
   mul(index)= v(20)/sig;
   co2_r_ex(index)=v(4);
   act_ex(index)=v(5);
   %prot_dmd(index)=sum(v(11:15));
   sur=sur+increment;
   x0 = v;
   end
   muary(1:loopcount,outloop)=mul;
   co2ary(1:loopcount,outloop)=co2_r_ex;
   actary(1:loopcount,outloop)=act_ex;
   %proary(1:loopcount,outloop)=prot_dmd;
end

% figure;
% h=zeros(1,5);
% colors =colormap(hsv(5));
% hold on
% for k=1:5
%     h(k)=plot(surl,muary(1:loopcount,k),'-','Color',colors(k,:),...
%         'DisplayName',sprintf('\\phi_z=%d%%',laczl(k)*100));
% end
% title({'Growth rates on varied substrate uptake bound','with different \\phi_z values'});
% xlabel('Substrate(Carbon) uptake rate bound(mmol OD600^{-1} h^{-1})');
% ylabel('Growth rate, \mu(h^{-1})');
% legend(h,'Location','NorthWest');
% hold off;
figure;
% h=zeros(3,5);
h=zeros(2,6);
colors =colormap(hsv(5));
hold on
for k=1:5
    h(1,k)=plot(muary(1:250,k),co2ary(1:250,k),'--','Color',colors(k,:),...
        'DisplayName',sprintf('co2:\\phi_z=%d%%',laczl(k)*100));
    h(2,k)=plot(muary(1:250,k),actary(1:250,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('act:\\phi_z=%d%%',laczl(k)*100));
%     h(3,k)=plot(muary(1:loopcount,k),proary(1:loopcount,k),'o-','Color',colors(k,:),...
%         'DisplayName',sprintf('pro:\\phi_z=%f',laczl(k)));    
end
h(1,6)=plot(MUL_NQ1389_phiZ_8,ACT_NQ1389_phiZ_8,'KS','DisplayName','act:\phi_z=8%');
h(2,6)=plot(MUL_NQ1389_phiZ_16,ACT_NQ1389_phiZ_16,'K^','DisplayName','act:\phi_z=16%');

% %parameters in basan's model:
% S_ac=1/3;
% %The paper says S_co2_r=1/6, which I believe it is wrong, and should be 1
% %instead.
% S_co2_r=1;
% %data from basan's model:
% mul=muary(1:10:loopcount,5);
% mu_co2=(1-phi0)*(1-lacZ/phi_max)/(sig/eps_f+b);
% mu_act=(1-phi0)*(1-lacZ/phi_max)/(sig/eps_r+b);
% co2_r_ex_d=(mul<=mu_co2).*(S_co2_r*(eps_r/e_r)*((1-phi0)*(1-lacZ/phi_max)-(sig/eps_f+b)*mul)/(1-(eps_r/eps_f)));
% act_ex_d=(mul>=mu_act).*(S_ac*(eps_f/e_f)*((sig/eps_r+b)*mul-(1-phi0)*(1-lacZ/phi_max))/((eps_f/eps_r)-1));
% h(1,7)=plot(mul,co2_r_ex_d./3,'o','DisplayName','Basan CO2:\phi_z=16%');
% h(2,7)=plot(mul,act_ex_d./3,'.','DisplayName','Basan act:\phi_z=16%');

legend(h(:),'Location','NorthWest');
title({'Growth-rate dependence of acetate excretion and respiratory CO_2 production','with different \phi_z values'});
xlabel('Growth rate, \mu(h^{-1})');
ylabel('Flux (mmol OD600^{-1} h^{-1})');
hold off;

%(2)energy dissipation
MUL_NQ1313_WT_LacY=[0.66548, 0.775801, 0.868327, 0.94306];
ACT_NQ1313_WT_LacY=[0., 0.596774, 1.48387, 1.95161];
MUL_NQ1313_Leaky_LacY=[0.491103, 0.562278, 0.629893, 0.661922, 0.672598];
ACT_NQ1313_Leaky_LacY=[0.741935, 1.82258, 2.20968, 2.66129, 2.90323];

beta=28.5;
alpha=1;
% eps_r=390;
% eps_f=750;
% e_r=4.4;
% e_f=2;
% b=0.13;
sig=45.7;
phi0=0.81;
phi_max=0.42;
lacZ=0;

eng_dis0=0;
kres=3.8;
kfer=41;
kribo=80;

loopcount=500;

eng_disl=zeros(1,5);
for outloop=1:5
   eng_dis=eng_dis0+(outloop-1)*16; 
   %lacZ=lacz0; 
   eng_disl(outloop)=eng_dis;
   mul=zeros(1,loopcount);
   %engdil=zeros(1,loopcount);
   act_ex=zeros(1,loopcount);
   co2_r_ex=zeros(1,loopcount);
   %prot_dmd=zeros(1,loopcount);
   %rp_ratio=zeros(1,loopcount);

   surl = zeros(1,loopcount);
   sur = 1;
   increment = 0.2;
   
%    kres=eps_r/(e_r*beta*alpha*(phi_max));
%    kfer=eps_f/(e_f*beta*alpha*(phi_max));
%    kribo=1/(b*(phi_max));
   x0 = double(ones(21,1));
   nonlcon = @(v)confun(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis);
   for index=1:loopcount
   surl(index)=sur;
   ub(1)=sur;
   %exitflag=0;
   %while exitflag ~= 1
   %    x0 = x0*1.01; 
       [v,fval,exitflag]=fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub,nonlcon,options);
   %end
   mul(index)= v(20)/sig;
   co2_r_ex(index)=v(4);
   act_ex(index)=v(5);
   %prot_dmd(index)=sum(v(11:15));
   sur=sur+increment;
   x0 = v;
   end
   muary(1:loopcount,outloop)=mul;
   co2ary(1:loopcount,outloop)=co2_r_ex;
   actary(1:loopcount,outloop)=act_ex;
   %proary(1:loopcount,outloop)=prot_dmd;
end

% figure;
% h=zeros(1,5);
% colors =colormap(hsv(5));
% hold on
% for k=1:5
%     h(k)=plot(surl,muary(1:loopcount,k),'-','Color',colors(k,:),...
%         'DisplayName',sprintf('energy dissipation=%d',eng_disl(k)));
% end
% title({'Growth rates on varied substrate uptake bound','with different energy dissipation values'});
% xlabel('Substrate(Carbon) uptake rate bound(mmol OD600^{-1} h^{-1})');
% ylabel('Growth rate, \mu(h^{-1})');
% legend(h,'Location','NorthWest');
% hold off;
figure;
% h=zeros(3,5);
%h=zeros(2,5);
h=zeros(2,6);
colors =colormap(hsv(5));
hold on
for k=1:5
    h(1,k)=plot(muary(1:200,k),co2ary(1:200,k),'--','Color',colors(k,:),...
        'DisplayName',sprintf('co2:energy dissipation=%d',eng_disl(k)));
    h(2,k)=plot(muary(1:200,k),actary(1:200,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('act:energy dissipation=%d',eng_disl(k)));
%     h(3,k)=plot(muary(1:loopcount,k),proary(1:loopcount,k),'o-','Color',colors(k,:),...
%         'DisplayName',sprintf('pro:\\phi_z=%f',laczl(k)));
end
h(1,6)=plot(MUL_NQ1313_WT_LacY,ACT_NQ1313_WT_LacY,'KS','DisplayName','act:LacY WT');
h(2,6)=plot(MUL_NQ1313_Leaky_LacY,ACT_NQ1313_Leaky_LacY,'K^','DisplayName','Leaky LacY');
legend(h(:),'Location','NorthWest');
title({'Growth-rate dependence of acetate excretion and respiratory CO_2 production','with different energy dissipation values'});
xlabel('Growth rate, \mu(h^{-1})');
ylabel('Flux (mmol OD600^{-1} h^{-1})');
hold off;

%(3)flh adn fli
%b
MUL_NQ1388_flhD=[1 0.97 0.9 0.82 0.77 0.6 1 0.97 0.98 0.96 0.96 0.94 1.02 1.01 0.98 0.93 0.89 0.86];
ACT_NQ1388_flhD=[1.53 0.32 0 0 0 0 2.82 2.29 2.18 2.62 1.44 1.31 2.47 2.31 0.85 0.22 -0.26 -0.09];

MUL_NQ1539_fliA=[1.01 0.93 0.85 0.81 0.65];
ACT_NQ1539_fliA=[2.34 0.23 -0.49 -0.33 0.05];
kres=5;
kfer=41;
kribo=80;
loopcount=500;
phi0_0=0.81;
phi0l=zeros(1,5);
for outloop=1:5
   phi0=phi0_0/(1.2^(outloop-1));
   phi0l(outloop)=phi0;
   mul=zeros(1,loopcount);
   act_ex=zeros(1,loopcount);
   surl = zeros(1,loopcount);
   sur = 1;
   increment = 0.2;

   nonlcon = @(v)confun(v,kres,kfer,kribo,alpha,beta,sig,phi0,phi_max,lacZ,eng_dis);
   for index=1:loopcount
   surl(index)=sur;
   ub(1)=sur;
   [v,fval]=fmincon(@objfun,x0,[],[],Aeq,beq,lb,ub,nonlcon,options);
   mul(index)= v(20)/sig;
   co2_r_ex(index)=v(4);
   act_ex(index)=v(5);
   sur=sur+increment;
   x0 = v;
   end
   muary(1:loopcount,outloop)=mul;
   actary(1:loopcount,outloop)=act_ex;
end

figure;
h=zeros(1,7);
colors =colormap(hsv(5));
hold on
for k=1:5
    h(1,k)=plot(muary(:,k),actary(:,k),'-','Color',colors(k,:),...
        'DisplayName',sprintf('act:\\phi_0=%f',phi0l(k)));
end
h(1,6)=plot(MUL_NQ1388_flhD,ACT_NQ1388_flhD,'KS','DisplayName','flhD');
h(1,7)=plot(MUL_NQ1539_fliA,ACT_NQ1539_fliA,'K^','DisplayName','fliA');

legend(h(:),'Location','NorthWest');
title({'Growth-rate dependence of acetate excretion and respiratory CO_2 production','with different \phi_0 values'});
xlabel('Growth rate, \mu(h^{-1})');
ylabel('Flux (mmol OD600^{-1} h^{-1})');
hold off;