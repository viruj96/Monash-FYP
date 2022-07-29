% A Mathematical Modelling Study of the Effects of Air Expansion Inside the
% Brain on the Intracranial Pressure

% Monash University Malaysia - Final Year Project
% Written by Viruj BALA SOUPRAMANIEN (27273652)
% Last modified: 10/9/19

% Change in ICP when change in altitude follows exponential and logarithmic
% functions

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

% Log function coefficients
X=[322.514 355.0392 394.8271];

for n=1:length(X)
    for m=1:length(VIA0)
        for l=1:length(PICr)
            A=VIA0(m)*PIA0; % From Boyles Law: P1V1=P2V2 (A=V_IA*P_IA)
            Patm=@(t)PIA0*(1-a*X(n).*t./(t+1)).^b; % Atmospheric pressure
            dPatm=@(t)-a*b*PIA0*(X(n)./(t+1)-X(n).*...
                t./(t+1).^2).*(1-a*X(n).*t./(t+1)).^(b-1); % Rate of change of atmospheric pressure
            V=@(t,P)A./(Patm(t)+P); % Intracranial air volume
            dPdt=@(t,P)(K*P./(R*((Patm(t)+P).^2+K*P*A))).*...
                ((PICr(l)-P).*(Patm(t)+P).^2-A*R*dPatm(t)); % Equation 8
            dV=@(t,P)-A*(dPatm(t)+dPdt(t,P))./((Patm(t)+P).^2);
            [t,P]=ode45(dPdt,[0 dt(n)],PICr(l)); % Solving equation 8
            icp{l,m,n}=P/133.322; % Intracranial pressure (mm Hg)
            alt{l,m,n}=X(n)*log(t+1)/0.3048; % Altitude (ft)
            iav{l,m,n}=V(t,P)*1e6; % Intracranial air volume (ml)
            delV{l,m,n}=dV(t,P)*1e9; % Rate of change of intracranial volume (ul/s)
            time{l,m,n}=t/60; % Time taken to reach 8000 ft (min)
            % Anderson Model
            V_Anderson=@(t)A./Patm(t); % Intracranial air volume
            dV_Anderson=@(t)-A*dPatm(t)./((Patm(t)).^2);
            dPdt_Anderson=@(t,P)(K*P/R)*(PICr(l)...
                +R*dV_Anderson(t)-P); % Equation 8
            [t,P]=ode45(dPdt_Anderson,[0 dt(n)],PICr(l)); % Solving equation 8
            P_IC{l,m,n}=P/133.322; % Intracranial pressure (mm Hg)
            H{l,m,n}=X(n)*log(t+1)/0.3048; % Altitude (ft)
            V_IA{l,m,n}=V_Anderson(t)*1e6; % Intracranial air volume (ml)
            rate{l,m,n}=dV_Anderson(t)*1e9; % Rate of change of intracranial volume (ul/s)
        end
    end
end

% Plots

[u,v,w]=size(icp);
PICr=PICr/133.322;
VIA0=VIA0/1e-6;
X=[2288 2824 3641];
lines={'r','b','g';'r--','b--','g--'};

for x=1:u
    for y=1:v
        if x==1
            c=1;
        else
            c=10;
        end
        for z=1:w
            figure(c)
            hold on
            plot(alt{x,y,z},iav{x,y,z},lines{1,y},'DisplayName',...
                ['Proposed Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            plot(H{x,y,z},V_IA{x,y,z},lines{2,y},'DisplayName',...
                ['Anderson Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            %title(sprintf('Change in Air Volume with Change in Altitude of %s*log(t+1) ft (P_{IC}_{_r} = %s mm Hg)',num2str(X(z)),num2str(PICr(x))))
            text(H{x,y,z}(end),V_IA{x,y,z}(end),...
                num2str(V_IA{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','bottom')
            text(H{x,y,z}(end),iav{x,y,z}(end),...
                num2str(iav{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','top')
            
            figure(c+1)
            hold on
            plot(alt{x,y,z},delV{x,y,z},lines{1,y},'DisplayName',...
                ['Proposed Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            plot(H{x,y,z},rate{x,y,z},lines{2,y},'DisplayName',...
                ['Anderson Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            %title(sprintf('Rate of Air Expansion against Change in Altitude of %s*log(t+1) ft (P_{IC}_{_r} = %s mm Hg)',num2str(X(z)),num2str(PICr(x))))
            
            figure(c+2)
            hold on
            plot(alt{x,y,z},icp{x,y,z},lines{1,y},'DisplayName',...
                ['Proposed Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            plot(H{x,y,z},P_IC{x,y,z},lines{2,y},'DisplayName',...
                ['Anderson Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            %title(sprintf('Change in ICP with Change in Altitude of %s*log(t+1) ft (P_{IC}_{_r} = %s mm Hg)',num2str(X(z)),num2str(PICr(x))))
            
            % Cleaning up presentation
            if z==1
                figure(c)
                axis([0 8500 0 70])
                if x==1
                    figure(c+2)
                    axis([0 13000 10 12.5])
                else
                    figure(c+2)
                    axis([0 12000 20 24.5])
                end
            elseif z==2
                figure(c)
                axis([0 8500 0 70])
                if x==1
                    figure(c+2)
                    axis([0 14000 10 12.5])
                else
                    figure(c+2)
                    axis([0 14000 20 25])
                end
            else
                figure(c)
                axis([0 8500 0 70])
                if x==1
                    figure(c+2)
                    axis([0 15000 10 13])
                else
                    figure(c+2)
                    axis([0 14000 20 25.5])
                end
            end
            
            figure(c)
            xlabel('Altitude [ft]','Fontweight','bold')
            ylabel('V_{IA} [ml]','Fontweight','bold')
            legend('-DynamicLegend','Location','best',...
                'Fontweight','bold')
            grid on
            grid minor
            figure(c+1)
            xlabel('Altitude [ft]','Fontweight','bold')
            ylabel('{dV_{IA}}/{dt} [\mul s^{-1}]','Fontweight','bold')
            legend('-DynamicLegend','Location','best',...
                'Fontweight','bold')
            grid on
            grid minor
            figure(c+2)
            xlabel('Altitude [ft]','Fontweight','bold')
            ylabel('Intracranial Pressure [mm Hg]','Fontweight','bold')
            legend('-DynamicLegend','Location','best',...
                'Fontweight','bold')
            grid on
            grid minor
            c=c+3;
        end
    end
end

% Exponential

PICr=PICr*133.322;
VIA0=VIA0*1e-6;

% Exp function coefficients
Y=[141.49 223.042 375.228];
Z=[0.0015 0.0026 0.0042];

for n=1:length(Y)
    for m=1:length(VIA0)
        for l=1:length(PICr)
            A=VIA0(m)*PIA0; % From Boyles Law: P1V1=P2V2 (A=V_IA*P_IA)
            Patm=@(t)PIA0*(1-a*Y(n)*Z(n).*...
                (exp(Z(n).*t-1)).*t).^b; % Atmospheric pressure
            dPatm=@(t)-a*b*PIA0*(Y(n)*Z(n).*...
                (exp(Z(n).*t-1))+Y(n)*Z(n).^2.*...
                (exp(Z(n).*t-1)).*t).*(1-a*Y(n)*...
                Z(n).*(exp(Z(n).*t-1)).*t).^(b-1); % Rate of change of atmospheric pressure
            V=@(t,P)A./(Patm(t)+P); % Intracranial air volume
            dPdt=@(t,P)(K*P./(R*((Patm(t)+P).^2+K*P*A))).*...
                ((PICr(l)-P).*(Patm(t)+P).^2-A*R*dPatm(t)); % Equation 8
            dV=@(t,P)-A*(dPatm(t)+dPdt(t,P))./((Patm(t)+P).^2);
            [t,P]=ode45(dPdt,[0 dt(n)],PICr(l)); % Solving equation 8
            icp{l,m,n}=P/133.322; % Intracranial pressure (mm Hg)
            alt{l,m,n}=Y(n).*(exp(Z(n).*t)-1)/0.3048; % Altitude (ft)
            iav{l,m,n}=V(t,P)*1e6; % Intracranial air volume (ml)
            delV{l,m,n}=dV(t,P)*1e9; % Rate of change of intracranial volume (ul/s)
            time{l,m,n}=t/60; % Time taken to reach 8000 ft (min)
            % Anderson Model
            V_Anderson=@(t)A./Patm(t); % Intracranial air volume
            dV_Anderson=@(t)-A*dPatm(t)./((Patm(t)).^2);
            dPdt_Anderson=@(t,P)(K*P/R)*(PICr(l)+...
                R*dV_Anderson(t)-P); % Equation 8
            [t,P]=ode45(dPdt_Anderson,[0 dt(n)],PICr(l)); % Solving equation 8
            P_IC{l,m,n}=P/133.322; % Intracranial pressure (mm Hg)
            H{l,m,n}=Y(n).*(exp(Z(n).*t)-1)/0.3048; % Altitude (ft)
            V_IA{l,m,n}=V_Anderson(t)*1e6; % Intracranial air volume (ml)
            rate{l,m,n}=dV_Anderson(t)*1e9; % Rate of change of intracranial volume (ul/s)
        end
    end
end

% Plots

PICr=PICr/133.322;
VIA0=VIA0/1e-6;
Y=[464 732 1231];
Z=[0.09 0.16 0.25];

for x=1:u
    for y=1:v
        if x==1
            c=19;
        else
            c=29;
        end
        for z=1:w
            figure(c)
            hold on
            plot(alt{x,y,z},iav{x,y,z},lines{1,y},'DisplayName',...
                ['Proposed Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            plot(H{x,y,z},V_IA{x,y,z},lines{2,y},'DisplayName',...
                ['Anderson Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            %title(sprintf('Change in Air Volume with Change in Altitude of %s*exp(%s*t-1) ft (P_{IC}_{_r} = %s mm Hg)',num2str(Y(z)),num2str(Z(z)),num2str(PICr(x))))
            text(H{x,y,z}(end),V_IA{x,y,z}(end),...
                num2str(V_IA{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','bottom')
            text(H{x,y,z}(end),iav{x,y,z}(end),...
                num2str(iav{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','top')
            
            figure(c+1)
            hold on
            plot(alt{x,y,z},delV{x,y,z},lines{1,y},'DisplayName',...
                ['Proposed Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            plot(H{x,y,z},rate{x,y,z},lines{2,y},'DisplayName',...
                ['Anderson Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            %title(sprintf('Rate of Air Expansion against Change in Altitude of %s*exp(%s*t+1) ft (P_{IC}_{_r} = %s mm Hg)',num2str(Y(z)),num2str(Z(z)),num2str(PICr(x))))
            text(H{x,y,z}(end),rate{x,y,z}(end),...
                num2str(rate{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','bottom')
            text(H{x,y,z}(end),delV{x,y,z}(end),...
                num2str(delV{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','bottom')
            
            figure(c+2)
            hold on
            plot(alt{x,y,z},icp{x,y,z},lines{1,y},'DisplayName',...
                ['Proposed Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            plot(H{x,y,z},P_IC{x,y,z},lines{2,y},'DisplayName',...
                ['Anderson Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            %title(sprintf('Change in ICP with Change in Altitude of %s*exp(%s*t+1) ft (P_{IC}_{_r} = %s mm Hg)',num2str(Y(z)),num2str(Z(z)),num2str(PICr(x))))
            text(H{x,y,z}(end),P_IC{x,y,z}(end),...
                num2str(P_IC{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','bottom')
            text(H{x,y,z}(end),icp{x,y,z}(end),...
                num2str(icp{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','top')
            
            % Cleaning up presentation
            if z==1
                figure(c)
                axis([0 9000 0 75])
            elseif z==2
                figure(c)
                axis([0 9000 0 75])
            else
                figure(c)
                axis([0 9000 0 75])
            end
            
            figure(c)
            xlabel('Altitude [ft]','Fontweight','bold')
            ylabel('V_{IA} [ml]','Fontweight','bold')
            legend('-DynamicLegend','Location','best',...
                'Fontweight','bold')
            grid on
            grid minor
            figure(c+1)
            xlabel('Altitude [ft]','Fontweight','bold')
            ylabel('{dV_{IA}}/{dt} [\mul s^{-1}]','Fontweight','bold')
            legend('-DynamicLegend','Location','best',...
                'Fontweight','bold')
            grid on
            grid minor
            figure(c+2)
            xlabel('Altitude [ft]','Fontweight','bold')
            ylabel('Intracranial Pressure [mm Hg]','Fontweight','bold')
            legend('-DynamicLegend','Location','best',...
                'Fontweight','bold')
            grid on
            grid minor
            c=c+3;
        end
    end
end
