function [qd,Ud]=divergence(rho,GJw,b,e,c,CL_alpha,Lambda)
Lambda
% COMPUTES THE TORSO-BENDING DIVERGENCE OF THE WING THROUGH
% THE RITZ-GALERKIN METHOD

syms y
n_th=2; %number of shape functions for torsion
n_w=2; %number of shape functions for bending
%% TORSION SHAPE FUNCTIONS

for k=1:n_th
    N_th(k) = sin((k*pi/(2*b)) * y); %row vector
end


%% BENDING SHAPE FUNCTIONS

for k=1:n_w
    N_w(k) = sin((k*pi/b) * y); %row vector
end


%% BUILD THE MATRICES FOR THE GENERALIZED EIGENVALUE PROBLEM

%Kthth
integrand_1=(diff(N_th))'*(GJw)*(diff(N_th));
A=int(integrand_1,[0 b]);

integrand_2=N_th'*e*c*CL_alpha*N_th;
B=int(integrand_2,[0 b]);

%Kthw
integrand_3=N_th'*e*c*CL_alpha*tan(Lambda)*diff(N_w);
C=int(integrand_3,[0 b]);

%Kwth
integrand_4=N_w'*c*CL_alpha*N_th;
D=int(integrand_4,[0 b]);

%Kww
integrand_5=N_w'*c*CL_alpha*tan(Lambda)*diff(N_w);
E=int(integrand_5,[0 b]);

%% GENERALIZED EIGENVALUE PROBLEM
A_tilde=[A , zeros(size(C,1),size(C,2));
        zeros(size(D,1),size(D,2)), zeros(size(E,1),size(E,2))];
B_tilde=[B,-C;
        -D, -E];
q=eig(double(A_tilde),double(B_tilde));
q(q<0)=nan;
U=sqrt(2*q./rho);
%% DIVERGENCE
qd=min(q);
Ud=sqrt(2*qd./rho);