% A Mathematical Modelling Study of the Effects of Air Expansion Inside the
% Brain on the Intracranial Pressure

% Monash University Malaysia - Final Year Project
% Written by Viruj BALA SOUPRAMANIEN (27273652)
% Last modified: 29/5/19

% Model 1
% Replicating results from: Andersson et al. (2003) - Air Transport of
% Patients with Intracranial Air: Computer Model of Pressure Effects

clc; clear; close all;

% Table 2 values
PICr=input('Initial resting pressure (mm Hg)? ')*133.322; % Resting pressure (Pa)
PVI=12.6*1e-6; % Pressure-volume index (m^3)
R=16.1*8.0124e9; % Outflow resistance (Pa/(m^3.s))
VIA0=[10 20 30 40 50]*1e-6; % Initial intracranial volume (m^3)
dHdt=input('Ascension rate (ft/min)? ')*0.3048/60; % Rate of ascension (m/s)
Hmax=8000*0.3048; % Maximum altitude (m)

% Numerical constants
a=2257e-8; % Alpha
b=5.264; % Beta
K=1/(0.4343*PVI); % Mathematical constant

for k=1:length(dHdt)
    dt(k)=Hmax/dHdt(k); % Time taken to reach 8000 ft
    for j=1:length(VIA0)
        V=@(t)VIA0(j)./(1-a.*dHdt(k).*t).^b; % Equation 13
        dVdt=@(t)a*b*VIA0(j).*dHdt(k)./((1-a.*dHdt(k).*t).^(b+1)); % Rate of change of intracranial air equation
        for i=1:length(PICr)
            dPdt=@(t,P)(K*P/R)*(PICr(i)+R*dVdt(t)-P); % Equation 8
            [t,P]=ode45(dPdt,[0 dt(k)],PICr(i)); % Solving equation 8
            H{i,j,k}=dHdt(k)*t/0.3048; % Cabin altitude (ft)
            P_IC{i,j,k}=P/133.322; % Intracranial pressure (mm Hg)
            V_IA{i,j,k}=V(t)*1e6; % Intracranial air volume (ml)
            rate{i,j,k}=dVdt(t)*1e9; % Rate of change of volume (ul/s)
        end
    end
end
for z=1:length(dHdt)
    for x=1:length(PICr)
        for y=1:length(VIA0)
            disp(P_IC{x,y,z}(end))
        end
    end
end
%% Plots
clear c d e u v w x y z L X Y lines figs; close all;


PICr=PICr/133.322;
VIA0=VIA0/1e-6;
dHdt=dHdt*60/0.3048;


[u,v,w]=size(P_IC);


for x=1:u
    for y=1:v
        for z=1:w
            
            figure(z)
            hold on
            plot(H{x,y,z},P_IC{x,y,z},'DisplayName',...
                [num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            title(sprintf('Change in ICP with Altitude at %s ft/min (P_{IC}_{_r} = %s mm Hg)',num2str(dHdt(z)),num2str(PICr(x))))
            legend('-DynamicLegend','Location','best')
            xlabel('Altitude [ft]')
            ylabel('Intracranial Pressure [mm Hg]')
            grid on
            grid minor
            
        end
    end
end
