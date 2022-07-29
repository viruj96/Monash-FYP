% A Mathematical Modelling Study of the Effects of Air Expansion 
% Inside the Brain on the Intracranial Pressure

% Monash University Malaysia - Final Year Project
% Written by Viruj BALA SOUPRAMANIEN (27273652)
% Last modified: 13/8/19

clc; clear; close all;

% Expressions for varying ascension rates

% Table 2 values
PICr=[10 20]*133.322; % Resting pressure (Pa)
VICr=[150 120]*1e-6; % Resting CSF volume
PVI=12.6*1e-6; % Pressure-volume index (m^3)
R=16.1*8.0124e9; % Outflow resistance (Pa/(m^3.s))
VIA0=[10 20 30]*1e-6; % Initial intracranial volume (m^3)
dHdt=[250 500 1000]*0.3048/60; % Rate of ascension (m/s)
Hmax=8000*0.3048; % Maximum altitude (m)

% Numerical constants
a=2257e-8; % Alpha
b=5.264; % Beta
K=1/(0.4343*PVI); % Mathematical constant
PIA0=101e3; % Absolute initial intracranial air (Pa) = P_atm at sea level

% Model 3
for n=1:length(dHdt)
    dt(n)=Hmax/dHdt(n); % Time taken to reach 8000 ft
    for m=1:length(VIA0)
        for l=1:length(PICr)
            A=VIA0(m)*PIA0; % From Boyles Law: P1V1=P2V2 (A=V_IA*P_IA)
            Patm=@(t)PIA0*(1-a*dHdt(n).*t).^b; % Atmospheric pressure
            dPatm=@(t)-a*b*PIA0*dHdt(n).*(1-a*dHdt(n).*t).^(b-1); % Rate of change of atmospheric pressure
            V=@(t,P)A./(Patm(t)+P); % Intracranial air volume
            dPdt=@(t,P)(K*P./(R*((Patm(t)+P).^2+K*P*A))).*((PICr(l)-P).*(Patm(t)+P).^2-A*R*dPatm(t)); % Equation 8
            dV=@(t,P)-A*(dPatm(t)+dPdt(t,P))./((Patm(t)+P).^2);
            [t,P]=ode45(dPdt,[0 dt(n)],PICr(l)); % Solving equation 8
            icp{l,m,n}=P/133.322; % Intracranial pressure (mm Hg)
            icv{l,m,n}=(log(P./PICr(l))/K+VICr(l))*1e6; % Intracranial volume (ml)
            iav{l,m,n}=V(t,P)*1e6; % Intracranial air volume (ml)
            iap{l,m,n}=(Patm(t)+P)/133.322; % Intracranial air pressure (mm Hg)
            delV{l,m,n}=dV(t,P)*1e9; % Rate of change of intracranial volume (ul/s)
            dicv{l,m,n}=icv{l,m,n}-VICr(l)*1e6; % Change in CSF volume (ml)
            dicp{l,m,n}=icp{l,m,n}-PICr(l)/133.322; % Change in intracranial pressure (mm Hg)
        end
    end
    time{n}=t/60; % Time taken to reach 8000 ft (min)
    h{n}=dHdt(n)*t/0.3048; % Altitude (ft)
end

% Time taken to reach 8000 ft @
T1=linspace(0,time{1}(end),41)'; % 250 ft/min
T2=linspace(0,time{2}(end),41)'; % 500 ft/min
T3=linspace(0,time{3}(end),41)'; % 1000 ft/min

% Logarithmic expression for ascension [H=A*log(t+1)]
Hlog1=@(t)2288*log(t+1); % 250 ft/min
Hlog2=@(t)2824*log(t+1); % 500 ft/min
Hlog3=@(t)3641*log(t+1); % 1000 ft/min

% Exponential expression for ascension [H=A*(exp(b*t)-1)]
Hexp1=@(t)464.2*(exp(0.0907*t)-1); % 250 ft/min
Hexp2=@(t)731.8*(exp(0.155*t)-1); % 500 ft/min
Hexp3=@(t)1231*(exp(0.2518*t)-1); % 1000 ft/min

% Plots
lines={'r','b','g';'r--','b--','g--';'rs-','bs-','gs-';'k','k--','ks-'};
figure
hold on
plot(T1,Hlog1(T1),'rs-',T1,Hexp1(T1),'r--'...
    ,T2,Hlog2(T2),'bs-',T2,Hexp2(T2),'b--'...
    ,T3,Hlog3(T3),'gs-',T3,Hexp3(T3),'g--')
for k=1:n
    plot(time{k},h{k},lines{1,k})
    Leg(k)=plot(NaN,NaN,lines{1,k});
    Leg(k+3)=plot(NaN,NaN,lines{4,k});
end
legend(Leg,'250 ft/min','500 ft/min','1000 ft/min',...
    'Linear','Exp','Log','Location','best','FontWeight','bold')
xlabel('Time (mins)','FontWeight','bold')
ylabel('Altitude (ft)','FontWeight','bold')
grid on
grid minor
