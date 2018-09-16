%%
%data on figure 1:
MUL_NCM3722_WT=[1.01 1.03 0.72 0.99 0.55 0.56 0.69 0.78 0.72 0.93 0.45 1.14];
ACT_NCM3722_WT=[2.32 2.84 0.07 2.67 -0.01 -0.01 0 -0.09 0 1.44 0 4.07];

MUL_NQ1243_UT=[0.58 0.64 0.82 0.95];
ACT_NQ1243_UT=[0 0.02 0.62 2.06];

MUL_NQ381_UT=[0.35 0.55 0.65 0.76 0.87 0.92];
ACT_NQ381_UT=[-0.01 -0.01 -0.01 -0.01 0.96 2.3];

MUL_NQ636_NQ638_NQ640_gm=[0.97 0.9 0.95];
ACT_NQ636_NQ638_NQ640_gm=[1.98 1.38 1.76];

MUL_NQ3722_7AA=[1.09 1.16 1.04 0.78 0.8 1.26];
ACT_NQ3722_7AA=[3.78 4.15 3.06 0.58 -0.01 4.18];

%Part 0: fit the data above:
beta=28.5;
alpha=1;
phi0=0.81;
phi_max=0.42;
lacZ=0;
sig=45.7;
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
%engdil=zeros(1,loopcount);
act_ex=zeros(1,loopcount);
co2_r_ex=zeros(1,loopcount);
%rp_ratio=zeros(1,loopcount);

surl = zeros(1,loopcount);
sur = 1;
increment = 0.2;
% 
kres=4.5;
kfer=80;
kribo=42;

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

plot(mul,act_ex,'r',...
     MUL_NCM3722_WT,ACT_NCM3722_WT,'KO',...
     MUL_NQ1243_UT,ACT_NQ1243_UT,'MO',...
     MUL_NQ381_UT,ACT_NQ381_UT,'M^',...
     MUL_NQ636_NQ638_NQ640_gm,ACT_NQ636_NQ638_NQ640_gm,'MS',...
     MUL_NQ3722_7AA,ACT_NQ3722_7AA,'KD');
title('Acetate excretion with varied growth rate (Batch Culture)');
xlabel('Growth rate, \mu(h^{-1})');
ylabel('Flux (mmol OD600^{-1} h^{-1})');
legend('ME:v_{EX\_acetate}','NCM3722 WT',...
       'NQ1243 Uptake Titration','NQ381 Uptake Titration',...
       'NQ636/NQ638/NQ640 glpK mutant','NQ3722 +7AA',...
       'Location','northwest');