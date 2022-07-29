% A Mathematical Modelling Study of the Effects of Air Expansion Inside the
% Brain on the Intracranial Pressure
% Effects of Temperature

% Monash University Malaysia - Final Year Project
% Written by Viruj BALA SOUPRAMANIEN (27273652)
% Last modified: 18/9/19

clc; clear all; close all;

% Parameters
PICr=[10 20]*133.322; % Resting pressure (Pa)
PVI=12.6*1e-6; % Pressure-volume index (m^3)
R=16.1*8.0124e9; % Outflow resistance (Pa/(m^3.s))
VIA0=[10 20 30]*1e-6; % Initial intracranial volume (m^3)
T0=294.15; % Initial body temperature (K)
Tf=310.15; % Final body temperature (K)
dT=1.55*34.5; % Rate of change of body temperature (K/s)

% Mathematical constant
K=1/(0.4343*PVI);

dt=(Tf-T0)/dT; % Time taken to warm body (s)

for m=1:length(VIA0)
    A=VIA0(m)./T0; % From Charles Law: V1/T1=V2/T2 (A=V_IA/T_f)
    dVdt(m)=A.*dT; % Rate of change of intracranial air volume (m^3/s)
    for l=1:length(PICr)
        dPdt=@(t,P)(K*P/R)*(PICr(l)+R*dVdt(m)-P); % Equation 8
        [t,P]=ode45(dPdt,[0 dt],PICr(l)); % Solving equation 8
        icp{l,m}=P/133.322; % Intracranial pressure (mm Hg)
        iav{l,m}=A.*dT.*t*1e6; % Intracranial air volume (ml)
    end
    delV(m)=dVdt(m)*1e9; % Rate of change of intracranial air volume (ul/s)
end
T=T0+dT*t-273.15; % Change in body temperature (K)

[v,w]=size(icp);
PICr=PICr/133.322;
VIA0=VIA0/1e-6;
lines={'r','b','g'};

figure(1)
plot(t,T,'k')
title('Change in Temperature with Time')
xlabel('Time [s]')
ylabel('Temperature [\circC]')
grid on
grid minor
text(0.1496,29.5,['T(t) = ',num2str(T0),' + ',num2str(dT),'*t [K]'],...
    'HorizontalAlignment','center','VerticalAlignment','bottom',...
    'Rotation',35,'FontWeight','bold')

for x=1:v
    for y=1:w
        if x==1
            c=2;
        else
            c=3;
        end
        figure(c)
        hold on
        plot(T,icp{x,y},lines{y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0}'])
        title(sprintf('Change in ICP with Change in Temperature (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
        xlabel('Temperature [\circC]')
        ylabel('Intracranial Pressure [mm Hg]')
        legend('-DynamicLegend','Location','best')
        grid on
        grid minor
        
        figure(c+2)
        hold on
        plot(T,iav{x,y},lines{y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0}'])
        title(sprintf('Change in Intracranial Air Volume with Temperature (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
        xlabel('Temperature [\circC]')
        ylabel('\DeltaV_{IA} [ml]')
        legend('-DynamicLegend','Location','best')
        grid on
        grid minor
    end
end