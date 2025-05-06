%% MATRICI DELLE EQUAZIONI DI MOTO NON LINEARI
% Questa funzione restituisce i termini destri delle equazioni di moto nella
% notazione state-space:
%       dtheta=dtheta
%       dphi=dphi
%       ddtheta=f(theta,phi,dtheta,dphi)
%       ddphi=f(theta,phi,dtheta,dphi)

% I valori delle costanti sono caricati dallo script parameters.mlx

function dq=stateEq_pend_inv(t,x,u,Tu)
parameters; % Carico i parametri

u = interp1(Tu,u,t);
theta = x(1);
phi = x(2);
dtheta = x(3);
dphi = x(4);
voltage = u;

%Mass Matrix
m11=(Jrod+Jenc+Jtau+mp*r^2+Jtheta_corr)+(mp*Lp^2)/3*sin(phi)^2;
m22=(mp*Lp^2)/3;
m12=-(mp*Lp*r)/2*cos(phi);
M=[m11, m12;m12, m22];

%Non linear terms
fNL=[(1/3)*mp*Lp^2*sin(2*phi)*dtheta*dphi+(1/2)*mp*Lp*r*sin(phi)*dphi^2;...
    -(1/6)*mp*Lp^2*sin(2*phi)*dtheta^2-(1/2)*g*mp*Lp*sin(phi)];
%Linear terms
fL=[Ctheta*dtheta+Kwire*theta;...
    Cenc*dphi];
%Relation Voltage-Torque of motor DC
tau=K_tau/R_tau*voltage; %L'effetto della back-EMF -K^2/R*dtheta è già contenuto in Ctheta
ftau=[tau;0];
%Vector with linear and non linear terms (also input action moved to left terms)
vet_fun=fNL+fL-ftau;

%State space rapresentation of system
dq=[dtheta;dphi;-M\vet_fun];
end