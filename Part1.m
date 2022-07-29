% A Mathematical Modelling Study of the Effects of Air Expansion Inside the
% Brain on the Intracranial Pressure

% Monash University Malaysia - Final Year Project
% Written by Viruj BALA SOUPRAMANIEN (27273652)
% Last modified: 29/5/19

% Anderson et al. Model
% Replicating results from: Andersson et al. (2003) - Air Transport of
% Patients with Intracranial Air: Computer Model of Pressure Effects

clc; clear; close all;

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
            V_IC{i,j,k}=(log(P./PICr(i))/K+VICr(i))*1e6; % Intracranial volume (ml)
            V_IA{i,j,k}=V(t)*1e6; % Intracranial air volume (ml)
            P_IA{i,j,k}=(101000*(1-a*dHdt(k).*t).^b)/133.322;
            rate{i,j,k}=dVdt(t)*1e9; % Rate of change of volume (ul/s)
        end
    end
end

% Proposed Model
% Verifying Anderson model using function of intracranial pressure [P_IC(t)]

clear A Patm dPatm l m n o t P dPdt V dV dicv icp icv iav iap delV alt;

PIA0=101e3; % Initial intracranial air pressure

for n=1:length(dHdt)
    for m=1:length(VIA0)
        for l=1:length(PICr)
            A=VIA0(m)*PIA0; % From Boyles law: P1V1=P2V2 (A=V_IA*P_IA)
            Patm=@(t)PIA0*(1-a*dHdt(n).*t).^b; % Atmospheric pressure
            dPatm=@(t)-a*b*PIA0*dHdt(n).*(1-a*dHdt(n).*t).^(b-1); % Rate of change of atmospheric pressure
            V=@(t,P)A./(Patm(t)+P);
            dPdt=@(t,P)(K*P./(R*((Patm(t)+P).^2+K*P*A))).*...
                ((PICr(l)-P).*(Patm(t)+P).^2-A*R*dPatm(t)); % Equation 8
            dV=@(t,P)-A*(dPatm(t)+dPdt(t,P))./((Patm(t)+P).^2);
            [t,P]=ode45(dPdt,[0 dt(n)],PICr(l)); % Solving equation 8
            icp{l,m,n}=P/133.322; % Intracranial pressure (mm Hg)
            icv{l,m,n}=(log(P./PICr(l))/K+VICr(l))*1e6; % Intracranial volume (ml)
            alt{l,m,n}=dHdt(n)*t/0.3048; % Altitude (ft)
            iav{l,m,n}=V(t,P)*1e6; % Intracranial air volume (ml)
            iap{l,m,n}=(Patm(t)+P)/133.322;
            delV{l,m,n}=dV(t,P)*1e9; % Rate of change of intracranial volume (ul/s)
        end
    end
end

% Plots

if PICr(1)~=10 && VIA0(1)~=10 && dHdt(1)~=250
    PICr=PICr/133.322;
    VIA0=VIA0/1e-6;
    dHdt=dHdt*60/0.3048;
end

[u,v,w]=size(icp);
lines={'r','b','g';'rs:','bs:','gs:';'r--','b--','g--'};

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
            plot(H{x,y,z},V_IA{x,y,z},lines{1,y},'DisplayName',...
                ['Anderson Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            plot(alt{x,y,z},iav{x,y,z},lines{3,y},'DisplayName',...
                ['Proposed Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            %title(sprintf('Change in Air Volume with Altitude at %d ft/min (P_{IC}_{_r} = %d mm Hg)',dHdt(z),PICr(x)))
            legend('-DynamicLegend','Location','best',...
                'FontWeight','bold')
            text(H{x,y,z}(end),V_IA{x,y,z}(end),...
                num2str(V_IA{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','bottom')
            text(H{x,y,z}(end),iav{x,y,z}(end),...
                num2str(iav{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','top')
            
            figure(c+1)
            hold on
            plot(H{x,y,z},rate{x,y,z},lines{1,y},'DisplayName',...
                ['Anderon Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            plot(alt{x,y,z},delV{x,y,z},lines{3,y},'DisplayName',...
                ['Proposed Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            %title(sprintf('Rate of Air Expansion against Altitude at %d ft/min (P_{IC}_{_r} = %d mm Hg)',dHdt(z),PICr(x)))
            legend('-DynamicLegend','Location','best',...
                'FontWeight','bold')
            text(H{x,y,z}(end),rate{x,y,z}(end),...
                num2str(rate{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','bottom')
            text(H{x,y,z}(end),delV{x,y,z}(end),...
                num2str(delV{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','top')
            
            figure(c+2)
            hold on
            plot(H{x,y,z},P_IC{x,y,z},lines{1,y},'DisplayName',...
                ['Anderson Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            plot(alt{x,y,z},icp{x,y,z},lines{3,y},'DisplayName',...
                ['Proposed Model: ',num2str(VIA0(y)),'ml V_{IA}_{_0}'])
            %title(sprintf('Change in ICP with Altitude at %d ft/min (P_{IC}_{_r} = %d mm Hg)',dHdt(z),PICr(x)))
            legend('-DynamicLegend','Location','best',...
                'FontWeight','bold')
            text(H{x,y,z}(end),P_IC{x,y,z}(end),...
                num2str(P_IC{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','bottom')
            text(H{x,y,z}(end),icp{x,y,z}(end),...
                num2str(icp{x,y,z}(end)),'HorizontalAlignment',...
                'left','VerticalAlignment','top')
            
            % Cleaning up presentation
            if z==1
                figure(c)
                axis([0 8500 0 75])
                figure(c+1)
                axis([0 8500 0 12])
                if x==1
                    figure(c+2)
                    axis([0 8500 10 21])
                else
                    figure(c+2)
                    axis([0 8500 20 31])
                end
            elseif z==2
                figure(c)
                axis([0 8500 0 75])
                figure(c+1)
                axis([0 8500 0 25])
                if x==1
                    figure(c+2)
                    axis([0 8500 10 26])
                else
                    figure(c+2)
                    axis([0 8500 20 40])
                end
            else
                figure(c)
                axis([0 8500 0 75])
                figure(c+1)
                axis([0 8500 0 46])
                if x==1
                    figure(c+2)
                    axis([0 8500 10 35])
                else
                    figure(c+2)
                    axis([0 8500 20 50])
                end
            end
            
            figure(c)
            xlabel('Altitude [ft]','FontWeight','bold')
            ylabel('V_{IA} [ml]','FontWeight','bold')
            grid on
            grid minor
            figure(c+1)
            xlabel('Altitude [ft]','FontWeight','bold')
            ylabel('{dV_{IA}}/{dt} [\mul s^{-1}]','FontWeight','bold')
            grid on
            grid minor
            figure(c+2)
            xlabel('Altitude [ft]','FontWeight','bold')
            ylabel('ICP [mm Hg]','FontWeight','bold')
            grid on
            grid minor
            c=c+3;
        end
    end
end