% A Mathematical Modelling Study of the Effects of Air Expansion Inside the
% Brain on the Intracranial Pressure
% Change in ICP when change in altitude follows exponential and logarithmic
% functions

% Monash University Malaysia - Final Year Project
% Written by Viruj BALA SOUPRAMANIEN (27273652)
% Last modified: 18/9/19

clc; clear all; close all;

% Logarithmic

% Table 2 values
PICr=[10 20]*133.322; % Resting pressure (Pa)
PVI=12.6*1e-6; % Pressure-volume index (m^3)
R=16.1*8.0124e9; % Outflow resistance (Pa/(m^3.s))
VIA0=[10 20 30]*1e-6; % Initial intracranial volume (m^3)
dt=[32 16 8]*60; % Time to reach maximum altitude (s)

% Numerical constants
a=2257e-8; % Alpha
b=5.264; % Beta
K=1/(0.4343*PVI); % Mathematical constant
PIA0=101e3; % Absolute initial intracranial air (Pa) = P_atm at sea level

coeff1=[322.514 355.0392 394.8271]; % Rate of ascension following log function (m/s)

for n=1:length(coeff1)
    for m=1:length(VIA0)
        for l=1:length(PICr)
            A=VIA0(m)*PIA0; % From Boyles Law: P1V1=P2V2 (A=V_IA*P_IA)
            Patm=@(t)PIA0*(1-a*coeff1(n).*t./(t+1)).^b; % Atmospheric pressure
            dPatm=@(t)-a*b*PIA0*(coeff1(n)./(t+1)-coeff1(n).*t./(t+1).^2)...
                .*(1-a*coeff1(n).*t./(t+1)).^(b-1); % Rate of change of atmospheric pressure
            V=@(t,P)A./(Patm(t)+P); % Intracranial air volume
            dPdt=@(t,P)(K*P./(R*((Patm(t)+P).^2+K*P*A))).*((PICr(l)-P).*(Patm(t)+P).^2-A*R*dPatm(t)); % Equation 8
            dV=@(t,P)-A*(dPatm(t)+dPdt(t,P))./((Patm(t)+P).^2);
            [t,P]=ode45(dPdt,[0 dt(n)],PICr(l)); % Solving equation 8
            icp{l,m,n}=P/133.322; % Intracranial pressure (mm Hg)
            alt{l,m,n}=(coeff1(n)*log(t+1))/0.3048; % Altitude (ft)
            iav{l,m,n}=V(t,P)*1e6; % Intracranial air volume (ml)
            delV{l,m,n}=dV(t,P)*1e9; % Rate of change of intracranial volume (ul/s)
            time{l,m,n}=t/60; % Time taken to reach 8000 ft (min)
            datm{l,m,n}=dPatm(t)*60/133.322; % Change in atmospheric pressure
            dp{l,m,n}=dPdt(t,P)*60/133.322;
        end
    end
end

[u,v,w]=size(icp);
PICr=PICr/133.322;
VIA0=VIA0/1e-6;
coeff1=[2288 2824 3641];
lines={'r','b','g';'r--','b--','g--';'r.','b.','g.'};

for x=1:u
    for y=1:v
        for z=1:w
            if x==1
                c=1;
            else
                c=2;
            end
            figure(c)
            hold on
            plot(alt{x,y,z},icp{x,y,z},lines{z,y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0} H=',num2str(coeff1(z)),'*log(t+1) ft'])
            title(sprintf('Change in ICP with Change in Altitude (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
            xlabel('Altitude [ft]')
            ylabel('Intracranial Pressure [mm Hg]')
            legend('-DynamicLegend','Location','best')
            grid on
            grid minor
            
            figure(c+2)
            hold on
            plot(time{x,y,z},icp{x,y,z},lines{z,y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0} H=',num2str(coeff1(z)),'*log(t+1) ft'])
            title(sprintf('Change in ICP with Time (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
            xlabel('Time [mins]')
            ylabel('Intracranial Pressure [mm Hg]')
            legend('-DynamicLegend','Location','best')
            grid on
            grid minor
            
            figure(c+4)
            hold on
            plot(alt{x,y,z},datm{x,y,z},lines{z,y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0} H=',num2str(coeff1(z)),'*log(t+1) ft'])
            title(sprintf('Change with P_{atm} with Altitude (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
            xlabel('Altitude [ft]')
            ylabel('Atmospheric Pressure [mm Hg]')
            legend('-DynamicLegend','Location','best')
            grid on
            grid minor
            
            figure(c+6)
            hold on
            plot(alt{x,y,z},dp{x,y,z},lines{z,y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0} H=',num2str(coeff1(z)),'*log(t+1) ft'])
            title(sprintf('Rate of Change in ICP with Altitude (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
            xlabel('Altitude [ft]')
            ylabel('dP_{IC}/dt [mm Hg/mins]')
            legend('-DynamicLegend','Location','best')
            grid on
            grid minor
        end
    end
end