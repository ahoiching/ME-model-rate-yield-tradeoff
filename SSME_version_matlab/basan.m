function output = basan( mu,eps_r,eps_f,b )
%BASAN Summary of this function goes here
e_f=2;
sig=45.7;
phi0=0.81;

%parameters in basan's model:
S_ac=1/3;
mu_act=(1-phi0)/(sig/eps_r+b);
act=(mu>=mu_act).*(S_ac*(eps_f/e_f)*((sig/eps_r+b)*mu-(1-phi0))/((eps_f/eps_r)-1));
output=act;

end


