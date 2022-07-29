% A Mathematical Modelling Study of the Effects of Air Expansion Inside the
% Brain on the Intracranial Pressure

% Monash University Malaysia - Final Year Project
% Written by Viruj BALA SOUPRAMANIEN (27273652)
% Last modified: 29/9/19

% Effects of temperature on intracranial system with pneumocephalus

clc; clearvars -except dy; close all;

% Parameters
PICr=[10 20]*133.322; % Resting intracranial pressure (Pa)
PVI=12.6*1e-6; % Pressure-volume index (m^3)
R=16.1*8.0124e9; % Outflow resistance (Pa/(m^3.s))
VIA0=[10 20 30]*1e-6; % Initial intracranial volume (m^3)
Ti=[18 21 24]+273.15; % Initial body temperature (K)
dt=[5 8 15]; % Time taken to reach final body temperature (s)

% Mathematical constant
K=1/(0.4343*PVI);

% Coefficients
grad=[3.5394 2.9728 2.4502;2.3215 1.953 1.58125;...
    1.250733333 1.052533333 0.854733333];
C=[Ti;Ti;Ti];

for n=1:length(Ti)
    for m=1:length(VIA0)
        Temp=@(t)grad(m,n)*t+C(m,n);
        dT=grad(m,n);
        CL(m,n)=VIA0(m)./Ti(n); % From Charles Law: V1/T1=V2/T2 (A=V_IA/T_f) (m^3/K)
        dVdt=CL(m,n).*dT; % Rate of change of intracranial air volume (m^3/s)
        for l=1:length(PICr)
            dPdt=@(t,P)(K*P/R)*(PICr(l)+R*dVdt-P); % Equation 8
            [t,P]=ode45(dPdt,[0 dt(m)],PICr(l)); % Solving equation 8
            icp{l,m,n}=P/133.322; % Intracranial pressure (mm Hg)
        end
        iav{m,n}=CL(m,n).*Temp(t)*1e6; % Intracranial air volume (ml)
        T{m,n}=Temp(t)-273.15; % Body temperature (C)
        time{m,n}=t; % Time taken (s)
        delV{m,n}=dVdt*1e9; % Rate of change of intracranial air volume (ul/s)
    end
end


[u,v,w]=size(icp);
PICr=PICr/133.322;
VIA0=VIA0/1e-6;
Ti=Ti-273.15;
lines={'r','b','g';'r--','b--','g--';'r.','b.','g.'};

for x=1:u
    for y=1:v
        for z=1:w
            if x==1
                if z==1
                    c=1;
                elseif z==2
                    c=3;
                else
                    c=5;
                end
                figure(c)
                hold on
                plot(T{y,z},iav{y,z},lines{1,y},'DisplayName',...
                    ['V_{IA}_{0} = ',num2str(VIA0(y)),'ml'])
                legend('-DynamicLegend','Location','best'...
                    ,'FontWeight','bold')
                xlabel('Temperature [\circC]','FontWeight','bold')
                ylabel('V_{IA} [ml]','FontWeight','bold')
                grid on
                grid minor
                text(T{y,z}(end),iav{y,z}(end),...
                    num2str(iav{y,z}(end)),'HorizontalAlignment',...
                    'left','VerticalAlignment','bottom')
            else
                if z==1
                    c=6;
                elseif z==2
                    c=7;
                else
                    c=8;
                end
            end
            figure(c+1)
            hold on
            plot(T{y,z},icp{x,y,z},lines{1,y},'DisplayName',...
                ['V_{IA}_{0} = ',num2str(VIA0(y)),'ml'])
            %title(sprintf('Change in ICP with Change in Temperature (P_{IC}_{_r} = %s mm Hg)',num2str(PICr(x))))
            text(T{y,z}(end),icp{x,y,z}(end),...
                num2str(icp{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','bottom')
            legend('-DynamicLegend','Location','best'...
                ,'FontWeight','bold')
            xlabel('Temperature [\circC]','FontWeight','bold')
            ylabel('Intracranial Pressure [mm Hg]','FontWeight','bold')
            grid on
            grid minor
            c=c+2;
        end
    end
end