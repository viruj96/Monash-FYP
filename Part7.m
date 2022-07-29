clc;clear all; close all;

T_V10=readtable('Data_viruj.xlsx','Sheet','V10','Range','A1:D32');
T_V20=readtable('Data_viruj.xlsx','Sheet','V20','Range','A1:D32');
T_V30=readtable('Data_viruj.xlsx','Sheet','V30','Range','A1:D32');

t=T_V10.t;
T{1,1}=T_V10.T_18;
T{2,1}=T_V20.T_18;
T{3,1}=T_V30.T_18;
T{1,2}=T_V10.T_21;
T{2,2}=T_V20.T_21;
T{3,2}=T_V30.T_21;
T{1,3}=T_V10.T_24;
T{2,3}=T_V20.T_24;
T{3,3}=T_V30.T_24;

[u,v]=size(T);
Ti=[18 21 24];
VIA0=[10 20 30];
lines={'r','b','g';'r--','b--','g--';'rs-','bs-','gs-';'k','k--','ks-'};

for i=1:u
    for j=1:v
        if i==1
            c=1;
        elseif i==2
            c=2;
        else
            c=3;
        end
        figure(c)
        hold on
        plot(t,T{i,j},lines{1,j},'DisplayName',['T_{i} = ',num2str(Ti(j)),'\circC'])
        xlabel('Time [s]','FontWeight','bold')
        ylabel('Temperature [\circC]','FontWeight','bold')
        legend('-DynamicLegend','Location','best','FontWeight','bold')
        grid on
        grid minor
    end
end

for k=1:c
    figure(k)
    yline(37,'m','DisplayName','37\circC');
end