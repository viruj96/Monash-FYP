% A Mathematical Modelling Study of the Effects of Air Expansion Inside the
% Brain on the Intracranial Pressure
% Effects of Temperature

% Monash University Malaysia - Final Year Project
% Written by Viruj BALA SOUPRAMANIEN (27273652)
% Last modified: 20/9/19

clc; clear all; close all;

% Parameters
PICr=[10 20]*133.322; % Resting intracranial pressure (Pa)
PVI=12.6*1e-6; % Pressure-volume index (m^3)
R=16.1*8.0124e9; % Outflow resistance (Pa/(m^3.s))
VIA0=[10 20 30]*1e-6; % Initial intracranial volume (m^3)
Ti=294.15; % Initial body temperature (K)
Tf=310.15; % Final body temperature (K)
Ta=311.15; % Warming temperature (K)
h=8; % Blood convective heat transfer coefficient (W/(m^2.K))
rhob=1050; % Blood density (kg/m^3)
cb=3800; % Blood specific heat capacity (J/(kg.K))
k=0.582; % Blood thermal conductivity (W/(m.K))
H0=470000; % Enthalpy of net chemical reaction of O2 and glucose (J/mol)
Hb=28000; % Energy used to release oxygen from hemoglobin (J/mol)
O2=8; % O2 concentration in blood (mol/m^3)
OEF=1.01; % O2 extraction fraction
rho=1100; % Tissue density (kg/m^3)
CBF=1.25*1e-5; % Blood flow (m^3/s)
delT=0.87; % Blood and tissue temperature difference (K)
Vb=0.0014; % Cranial volume (m^3)
mb=1.4; % Mass of brain (kg)
Eg=Vb*rho*CBF*((H0-Hb)*O2*OEF-rhob*cb*delT)/mb; % Metabolic heat generation (W)

% Mathematical constant
K=1/(0.4343*PVI);

for m=1:length(VIA0)
    r(m)=(0.75*VIA0(m)/pi).^(1/3); % Radius of air bubble (m)
    % Applicability of lump capacitance model
    Bi=h*r/3/k; % Biot number
    if Bi<0.1
        As=4*pi*r(m)^2; % Surface area of air bubble (m^2)
        tcons=rhob*cb*VIA0(m)./(As*h); % Thermal time constant (s)
        a=1/tcons;
        b=(As*h*(Tf-Ti)+Eg)/(rhob*VIA0(m)*cb);
        dt(m)=-tcons*log((Tf-Ta-tcons*b)/(Ti-Ta-tcons*b)); % Time to warm body (s)
        Temp=@(t)Ta+(Ti-Ta)*exp(-a*t)+tcons*b*(1-exp(-a*t)); % Temperature distribution (K)
        dT=@(t)-a*(Ti-Ta)*exp(-a*t)+b*exp(-a*t); % Rate of change of temperature (K/s)
        A=VIA0(m)./Ti; % From Charles Law: V1/T1=V2/T2 (A=V_IA/T_f) (m^3/K)
        dVdt=@(t)A.*dT(t); % Rate of change of intracranial air volume (m^3/s)
        for l=1:length(PICr)
            dPdt=@(t,P)(K*P/R)*(PICr(l)+R*dVdt(t)-P); % Equation 8
            [t,P]=ode45(dPdt,[0 dt(m)],PICr(l)); % Solving equation 8
            icp{l,m}=P/133.322; % Intracranial pressure (mm Hg)
            iav{l,m}=A.*Temp(t)*1e6; % Intracranial air volume (ml)
            T{l,m}=Temp(t)-273.15; % Body temperature (C)
            time{l,m}=t/60; % Time taken (mins)
            delV{l,m}=dVdt(t)*1e9; % Rate of change of intracranial air volume (ul/s)
            DelT{l,m}=dT(t); % Rate of change of temperature with time (C/s)
        end
    else
        disp('Biot number>0.1 -> Lump capacitance model n/a')
    end
end

[v,w]=size(icp);
PICr=PICr/133.322;
VIA0=VIA0/1e-6;
lines={'r','b','g'};

for x=1:v
    for y=1:w
        figure(1)
        hold on
        plot(time{x,y},T{x,y},lines{y})
        if x==1
            c=2;
        else
            c=3;
        end
        figure(c)
        hold on
        plot(T{x,y},icp{x,y},lines{y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0}'])
        title(sprintf('Change in ICP with Change in Temperature (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
        xlabel('Temperature [\circC]')
        ylabel('Intracranial Pressure [mm Hg]')
        legend('-DynamicLegend','Location','best')
        grid on
        grid minor
        
        figure(c+2)
        hold on
        plot(T{x,y},iav{x,y},lines{y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0}'])
        title(sprintf('Change in Intracranial Air Volume with Temperature (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
        xlabel('Temperature [\circC]')
        ylabel('V_{IA} [ml]')
        legend('-DynamicLegend','Location','best')
        grid on
        grid minor
        
        figure(c+4)
        hold on
        plot(T{x,y},delV{x,y},lines{y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0}'])
        title(sprintf('Rate of Change of Intracranial Air Volume with Temperature (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
        xlabel('Temperature [\circC]')
        ylabel('dV_{IA}/dt [\mul/s]')
        legend('-DynamicLegend','Location','best')
        grid on
        grid minor
        
        figure(c+6)
        hold on
        plot(DelT{x,y},delV{x,y},lines{y},'DisplayName',[num2str(VIA0(y)),'ml V_{IA}_{_0}'])
        title(sprintf('Rate of Change of Temperature with Change in Volume (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
        xlabel('dT_{IA}/dt [\circC/mins]')
        ylabel('dV_{IA}/dt [\mul/s]')
        legend('-DynamicLegend','Location','best')
        grid on
        grid minor
    end
end

figure(1)
title(sprintf('Change in Temperature with Time'))
xlabel('Time [mins]')
ylabel('Temperature [\circC]')
legend({'10ml V_{IA}_{_0}','20ml V_{IA}_{_0}','30ml V_{IA}_{_0}'},'Location','best')
grid on
grid minor
