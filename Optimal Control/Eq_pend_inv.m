function dx=Eq_pend_inv(x,u)
parameters;

%Mass Matrix
m11 = (Jrod+Jenc+Jtau+mp*r^2+Jtheta_corr)+(mp*Lp^2)/3*sin(x(2))^2;
m22 = (mp*Lp^2)/3;
m12 = - (mp*Lp*r)/2*cos(x(2));
M=[m11, m12;m12, m22];

%Non linear terms
fNL=[(1/3)*mp*Lp^2*sin(2*x(2))*x(3)*x(4)+(1/2)*mp*Lp*r*sin(x(2))*x(4)^2;...
    -(1/6)*mp*Lp^2*sin(2*x(2))*x(3)^2-(1/2)*g*mp*Lp*sin(x(2))];

%Linear terms
fL=[Ctheta*x(3)+Kwire*x(1);...
    Cenc*x(4)];

%Relation Voltage-Torque of motor DC
tau=K_tau/R_tau*u; %L'effetto della back-EMF -K^2/R*dtheta è già contenuto in Ctheta
ftau=[tau;0];

%Vector with linear and non linear terms (also input action moved to left terms)
vet_fun=fNL+fL-ftau;

%State space rapresentation of system
dx= [x(3);x(4);-M\vet_fun];
end