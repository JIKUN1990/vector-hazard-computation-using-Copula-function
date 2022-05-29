clc
clear
% scenario eartqhuake information (the example case is the same in Baker et al(2021)Seismic hazard and risk analysis: P256)
M=6.5;
R=10;
Fault_Type=1; % strike-slip
Vs=500; % Vs30
arb=0;
lamda=0.01; % annual rate for (M,R).
%% Scalar PSHA for Sa(1.0s)
T=1.0 % studied period
[sa, sigma] = BJF_1997_horiz(M, R, T, Fault_Type, Vs, arb); % The GMPE propopsed by Boore et al(1997), which is coded by Jack Baker, 2/1/05
x=0.01:0.01:0.9;
for i=1:length(x)
    y(i)=normpdf((log(x(i))-log(sa))/sigma)./(x(i)*sigma); % mean rate density of Sa(1.0s)
end

for i=1:length(x)-1
    ARE(i)=trapz(x(i:end),y(i:end))*lamda; % Annual rate exceedance
end
Hazardcurve=[x(1:end-1)' ARE'];
plot(x(1:end-1)' , ARE') % hazard curve for Sa(1.0s)
%% Direct Vector PSHA (Sa(1.0s):0.42g; Sa(2.0s):0.20g)
T1=1.0;
[sa1, sigma1] = BJF_1997_horiz(M, R, T1, Fault_Type, Vs, arb);
T2=2.0;
[corr_req] = baker_jayaram_correlation(T1,T2);
[sa2, sigma2] = BJF_1997_horiz(M, R, T2, Fault_Type, Vs, arb);

x=0.42:0.01:1.0; % Sa(1.0s)=0.42g
for ii=1:length(x)
    y= 0.01:0.01:1;
    conditional_mean =log(sa2)+corr_req*sigma2/sigma1*(log(x(ii))-log(sa1));
    conditional_std  =sigma1*sqrt(1-corr_req^2);
    for i=1:length(y)
        yy(i)=normpdf((log(y(i))-conditional_mean)/conditional_std)./(y(i)*conditional_std) * normpdf((log(x(ii))-log(sa1))/sigma1)./(x(ii)*sigma1);
    end
    
    for i=1:length(y)-1
        Z(i)=trapz(y(i:end),yy(i:end))*lamda;
    end
    temp1(ii)=Z(20); % Sa(2.0s)=0.20g
end
MERs_direct_VPSHA=trapz(x,temp1);
%% MER calculation using Copula method (Sa(1.0s):0.42g; Sa(2.0s):0.20g)
T1=1.0 ;
[sa1, sigma1] = BJF_1997_horiz(M, R, T1, Fault_Type, Vs, arb);
P1= normcdf(log(0.42),log(sa1),sigma1); % P(Sa(1.0s)<0.42g|(M,R))

T2=2.0;
[corr_req] = baker_jayaram_correlation(T2,T1);
[sa2, sigma2] = BJF_1997_horiz(M, R, T2, Fault_Type, Vs, arb);
P2= normcdf((log(0.20)-log(sa2))/sigma2); % P(Sa(2.0s)<0.30g|(M,R))
% MERs("AND case");
z(1,1) = norminv(1-P1);
z(2,1) = norminv(1-P2);
MERs_copula = mvncdf(z,[0;0],[1 corr_req;corr_req 1])*lamda;

% MERa("OR case")
z(1,1) = norminv(P1);
z(2,1) = norminv(P2);
MERa_copula = 1-mvncdf(z,[0;0],[1 corr_req;corr_req 1])*lamda;

% MERk("KEN case")
z(1,1) = norminv(P1);
z(2,1) = norminv(P2);
KTsize=100;
N_IMs = 2;
t =mvncdf(z,[0;0],[1 corr_req;corr_req 1]);%defination of t
mm=1;

for ii=1:50
    
    sample=mvnrnd([0;0],[1 corr_req;corr_req 1],KTsize);
    
    k=0;
    v=[];
    for i=1:length(sample(:,1))
        v=mvncdf(sample(i,:)',[0;0],[1 corr_req;corr_req 1]);
        if t>=0 && v<=t
            k=k+1;
        else
        end
    end
    Kt=k/length(sample(:,1));
    if Kt==0
        Kt=nan;
    else
    end
    result(mm,2)=Kt;
    mm=mm+1;
    prob_unif(ii) = 1-Kt;
end

MERk_copula=[mean(prob_unif*lamda), std(prob_unif*lamda)]


