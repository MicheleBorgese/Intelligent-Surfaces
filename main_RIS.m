clear all
close all

% Computes the reflection coefficient of a RIS comprising an array of partches loased with varactor diodes.
%
% This function was developed as a part of the paper:
%
% Filippo Costa, Michele Borgese, “Electromagnetic Model of Reflective Intelligent Surfaces,” submitted to iEEE Trasactions On Wireless Communications.
%
% This is version 1.0 (Last edited: 2021-02-15)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.

% Frequency (GHz)
freq=(1:0.1:15); 

%Incidence angle (deg)
theta=0;

%Varactor capacitance range
Cap_vect=(0.1:0.1:0.5) ;
%Varactor inductance
Lvar=0.5e-9;

%RIS parameters
Dx=5e-3; %patch array periodicity along x-direction
Dy=Dx; %patch array periodicity along y-direction
wx=0.5e-3; %patch array gap width along x-direction
wy=wx; %patch array gap width along y-direction
sigma_copper=58.7*1e6; %copper conductivity
mu0 = 4 * pi * 1e-7;  %vacuum permeability 
delta=sqrt(1./(pi*freq*1e9*sigma_copper*mu0)); %skin depth
Rs=1./(sigma_copper*delta); %surface impedance
er1=4.4-1i*0.088; %substrate dielectric permittivity
d=1.2e-3; %substrate thickness (m)

Ncap = length(Cap_vect) ;

gamma_TE = zeros(length(freq), Ncap) ;
gamma_TM = zeros(length(freq), Ncap) ;

% computation of reflection coefficient
for kk=1:Ncap
[gamma_TE(:,kk),gamma_TM(:,kk)] = RIS_reflection(freq*1e9,Dx,wx,Dy,wy,Rs,d,er1,theta,Cap_vect(kk),Lvar);
end

%plot reflection coefficient amplitude and phase as a function of the varactor capacitance

index_plot=[1,2,3,4,5];

for pp=1:length(index_plot)
    
    legend_info{pp}=['Model - C = ' num2str(Cap_vect(index_plot(pp))) ' pF'];
    
    figure(1)
    plot(freq,20*log10(abs(gamma_TE(:,index_plot(pp)))),'linewidth',2); hold on
    
    figure(2)
    plot(freq,20*log10(abs(gamma_TM(:,index_plot(pp)))),'linewidth',2); hold on
    
    figure(3)
    plot(freq,180/pi*(angle(gamma_TE(:,index_plot(pp)))),'linewidth',2); hold on
    
    figure(4)
    plot(freq,180/pi*(angle(gamma_TM(:,index_plot(pp)))),'linewidth',2); hold on
    
end

min_reflection=-3;
font_size=18;

figure(1);
xlabel('Frequency (GHz)','interpreter','latex','Fontsize',font_size,'FontName','Times')
ylabel('$|\Gamma|$ (dB)','interpreter','latex','FontName','Times')
legend(legend_info,'interpreter','latex','FontName','Times','Fontsize',14,'Location','Southeast')
set(gca,'Fontsize',font_size)
axis([min(freq) max(freq) min_reflection 0])


figure(2)
xlabel('Frequency (GHz)','interpreter','latex','Fontsize',font_size,'FontName','Times')
ylabel('$|\Gamma|$ (dB)','interpreter','latex','FontName','Times')
legend(legend_info,'interpreter','latex','FontName','Times','Fontsize',14,'Location','Southeast')
set(gca,'Fontsize',font_size)
axis([min(freq) max(freq) min_reflection 0])


figure(3)
xlabel('Frequency (GHz)','interpreter','latex','Fontsize',font_size,'FontName','Times')
ylabel('$\angle{\Gamma}$ (deg)','interpreter','latex','FontName','Times')
legend(legend_info,'interpreter','latex','Fontsize',14,'FontName','Times','Location','Northeast')
set(gca,'Fontsize',font_size)
axis([min(freq) max(freq) -180 180])



figure(4)
xlabel('Frequency (GHz)','interpreter','latex','Fontsize',font_size,'FontName','Times')
ylabel('$\angle{\Gamma}$ (deg)','interpreter','latex','FontName','Times')
legend(legend_info,'interpreter','latex','Fontsize',14,'FontName','Times','Location','Northeast')
set(gca,'Fontsize',font_size)
axis([min(freq) max(freq) -180 180])




