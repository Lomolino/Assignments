clc
clear all 
close all

%% GENERAL DATA
% 0453
% DA=30
% CB=54

x1=1.2;
x2=0.6;

rho=1.225;
%% GEOMETRY DATA
b=6.1; %semi-wing span
c=3.05; %wing chord
Lambda=deg2rad(-30); %sweep angle

Lf=9.15; %fuselage length

cc=3.05; %canard chord
bc=1.525; %canard semi-wing span

%wing attached at (0,0) on the fuselage
%canard attached at (-4.525,o) on the fuselage

ca=1.4*c; %ailerons chord
ba=1/2*b; %ailerons semi-span

% the elastic axis of the wing is positioned at half of the chord
e=c/2;

%% STRUCTURAL PROPERTIES WING
EIw=x1*(4.5*10^6); 
GJw=x2*(7.0*10^6);

%% STRUCTURAL PROPERTIES FUSELAGE
GJf=12*10^6;

%% AERODYNAMICS
CL_alpha=2*pi; %both wing and canard
%aerodynamic center AC at 1/4 of the chord
%AERODYNAMICS OF THE AILERONS FROM SET 2 OF SLIDES

%% POINT 1
Lambda=linspace(deg2rad(-40),deg2rad(40),10);
qd=nan(length(Lambda));
Ud=nan(length(Lambda));
for i=1:length(Lambda)
    [qd(i),Ud(i)]=divergence(rho,GJw,b,e,c,CL_alpha,Lambda(i));
end

figure()
plot(rad2deg(Lambda),Ud);
grid minor
xlabel('Sweep Angle')
ylabel('Divergence Speed')
