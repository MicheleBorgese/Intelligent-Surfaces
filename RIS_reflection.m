function [gamma_TE,gamma_TM] = RIS_reflection(freq,Dx,wx,Dy,wy,Rs,d,er1,theta_deg,cap,Lvar)

%This function computes the reflection coefficient of a RIS comprising an array of patches loaded with varactor diodes.
%
% This function was developed as a part of the paper:
%
% Filippo Costa, Michele Borgese, “Electromagnetic Model of Reflective Intelligent Surfaces”, submitted to IEEE Transactions on Wireless Communications.
%
% This is version 1.0 (Last edited: 2021-02-15)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.
%
% INPUT:
% freq = analyzed frequency or frequency range
% theta_deg= angle of incidence (deg)
% cap= capacitance of the varactor (pF)
% Lvar= inductance of the varactor (H)
% Dx = periodicity of the periodic surface along x-direction (m)
% Dy = periodicity of the periodic surface along y-direction (m)
% wx = patch array gap width along x-direction (m)
% wy = patch array gap width along y-direction (m)
% Rs = surface impedance of the conductive material used to fabricate the patch array (ohm/sq)
% d=thcikness of the dielectric substrate (m)
% er1=dielectric permittivity of the dielectric substrate

% OUTPUT:
% gamma_TE, gamma_TM = reflection coefficient (complex) of the RIS for TE and TM polarization


% varactor capacitance
Cvar=cap*1e-12;

%effective permettivity
ereff=(er1+1)/2;

%incidence angle
theta=deg2rad(theta_deg);

% vacuum permittivity and permeability 
mu0 = 4 * pi * 1e-7;
eps0 = 8.85 * 1e-12;
c0=1/sqrt(eps0*mu0);

%wavelength
lambda=c0./freq;
omega=2*pi*freq;

%propagation constants
k0=omega*sqrt(eps0*mu0);
keff=k0*sqrt(ereff);
kt=k0*sin(theta); %transverse component
kz0=sqrt(k0.^2-kt.^2); %normal component in vacuum
kz1=sqrt(er1*k0.^2-kt.^2); %normal component in the substrate

%impedances
z0=sqrt(mu0/eps0);
Z0TE=omega.*mu0./kz0; %vacuum - TE
Z0TM=kz0./(omega*eps0); %vacuum - TM 
Z1TE=omega.*mu0./kz1; %substrate - TE - relation (4)
Z1TM=kz1./(omega*eps0*er1); %substrate - TM - relation (4)

%ohmic resistance
Rx=Rs*(Dx/(Dx-wx))^2;
Ry=Rs*(Dy/(Dy-wy))^2;

%patch capacitance
CTE_patch=2*Dy*eps0*ereff/pi*log(1/sin(pi*wy/(2*Dy)))*(1-k0.^2./keff.^2*1/2*sin(theta)^2); %relation (6)
CTM_patch=2*Dx*eps0*ereff/pi*log(1/sin(pi*wx/(2*Dx))); %relation (7)

%ground-patch capacitance
Cap_correction_x=2*eps0*Dx/pi.*log(1-exp(-4*pi*d./Dx));
Cap_correction_y=2*eps0*Dy/pi.*log(1-exp(-4*pi*d./Dx));

%patch admittance
Ypatch_TE=1i*omega.*(CTE_patch-Cap_correction_y);
Ypatch_TM=1i*omega.*(CTM_patch-Cap_correction_x);

%varactor impedance
Zvar=1i.*omega.*Lvar+1./(1i*omega.*Cvar); %relation (10)
Yvar=1./Zvar;

%Zsurf impedance
ZsurfTE=1./(Ypatch_TE+Yvar); 
ZsurfTM=1./(Ypatch_TM+Yvar);

%grounded substrate input impedance
Zd_TE = 1i*Z1TE.*tan(kz1*d) ; %relation (3)
Zd_TM = 1i*Z1TM.*tan(kz1*d) ; %relation (3)

%total input impedance
Zv_TE = 1./(1./ZsurfTE+1./Zd_TE); %relation (1)
Zv_TM = 1./(1./ZsurfTM+1./Zd_TM); %relation (1)

%reflection coefficient
gamma_TE=(Zv_TE-Z0TE)./(Zv_TE+Z0TE);
gamma_TM=(Zv_TM-Z0TM)./(Zv_TM+Z0TM);

end

