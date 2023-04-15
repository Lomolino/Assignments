close all
clear
clc

%% ASSIGNMENT 

% ABCD last four letter of the person code 10650070

A = 0;
B = 0;
C = 7;
D = 0;

DA = str2double(sprintf('%d%d', D, A)); 
CB = str2double(sprintf('%d%d', C, B));

if DA<50
    x1 = 1.2;
else 
    x1 = 0.7;
end

if CB<50
    x2 = 1.1;
else
    x2 = 0.6;
end

%% Dimnesions

b = 6.1;        % Wing semi-span
c = 3.05;       % Wing chord
lambda = -30*pi/180;    % Sweep angle
Lf = 9.15;      % Fuselage length
c_c = 3.05;     % Canard chord
b_c = 1.525;    % Canard semi-span

y_m = b/cos(lambda);
%% Wing-Canard Position
x_w = 0;
y_w = 0;

x_c = -4.525;
y_c = 0;

%% Ailerons

E = 0.25;
c_a = E*c;    % Ailerons chord
b_a = 0.5*b;    % Ailerons semi-span
y_a = b_a/cos(lambda);
e = 0.25*c;

beta = 1*pi/180;

%% Properties

Ea = 0.5;       % Elastic axis % position

EIw = x1*4.5*10^6;  % Bending stiffness
GJw = x2*7.0*10^6;  % Torsional stiffness
GJf = 12.0*10^6;    % Fuselage Torsional stiffness


Cl_alpha = 2*pi;     % Both for the wing and the canard

Cl_beta = 1/pi*(acos(1-2*E)+2*sqrt(E*(1-E)))*Cl_alpha;   % Aileron Cl
Cm_beta = -1/pi*(1-E)*sqrt(E*(1-E))*Cl_alpha;            % Aileron Cm


%% Shape Functions

dim_w = 5;
dim_t = 4;

y_w = linspace(0,y_m,dim_w);
y_t = linspace(0,y_m,dim_w);

for n = 2:(dim_w+1)
    Nw(:,n-1)=(y_w/b).^n;
    dNw(:,n-1)=n*(y_w/b).^(n-1);
    d2Nw(:,n-1)=n*(n-1)*(y_w/b).^(n-2);
end

for n = 1:(dim_t)
    Nt(:,n) = sin(n*pi/2*y_t/b);
    dNt(:,n) = n*pi/2/b*cos(n*pi/2*y_t/b);
end


%% Matrix creation

Kww_s = d2Nw'*EIw*d2Nw*y_m;
Ktt_s = dNt'*GJw*dNt*y_m;

Kww_a = ((-Nw'*c*Cl_alpha*sin(lambda)*cos(lambda)*dNw)+(dNw'*c*Cl_alpha*sin(lambda)^2*cos(lambda)*dNw))*y_m;
Kwt_a = ((Nw'*c*Cl_alpha*(cos(lambda))^2*Nt)-(dNw'*c*Cl_alpha*(cos(lambda))^2*sin(lambda)*Nt))*y_m;
Ktw_a = -Nt'*e*c*Cl_alpha*sin(lambda)*(cos(lambda))^2*dNw*y_m;
Ktt_a = Nt'*e*c*Cl_alpha*cos(lambda)^3*Nt*y_m;

Ks = [Ktt_s zeros(dim_t,dim_w); zeros(dim_w,dim_t) Kww_s];
Ka = [Ktt_a Ktw_a; Kwt_a Kww_a];

Fw = (Nw'*c*Cl_beta*cos(lambda)-dNw'*e*c*Cl_beta*sin(lambda)*cos(lambda)-(dNw'*c^2*Cm_beta*cos(lambda)*sin(lambda)))*(y_m-y_a);
Ft = ((Nt'*e*c*Cl_beta*cos(lambda)^2)+(Nt'*c^2*Cm_beta*cos(lambda)^2))*(y_m-y_a);



%% Divergence Dynamic pressure 

qd = eig(Ks, Ka);
qd = qd(qd>=0);

v = sqrt(abs(qd)*2/1.225);

[v_d, i] = min(v)
q_d = qd(i)

%figure(1)
%plot(v)

q = 1;

K = (Ks-q.*Ka);
F = [Ft;Fw];

z_t = inv(K)*F;
z = z_t*beta*q;


L = c*Cl_alpha*(Nt*z_t'*cos(lambda)-dNw*z_t'*sin(lambda))*cos(lambda)*y_m+c*Cl_beta*cos(lambda)*(y_m-y_a);



