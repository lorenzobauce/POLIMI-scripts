function fu_resub = fu_jacobian(in1,u,Jrod,Jenc,Jtau,mp,r,Jtheta_corr,Lp,g,K_tau,R_tau,Ctheta,Kwire,Cenc)
parameters;
x2 = in1(2,:);
t2 = sin(x2);
t3 = Lp.^2;
t4 = r.^2;
t6 = Jenc.*1.2e+1;
t7 = Jrod.*1.2e+1;
t8 = Jtau.*1.2e+1;
t9 = Jtheta_corr.*1.2e+1;
t10 = 1.0./R_tau;
t5 = t2.^2;
t11 = mp.*t4.*3.0;
t12 = mp.*t3.*t5.*4.0;
t13 = mp.*t4.*t5.*9.0;
t14 = t6+t7+t8+t9+t11+t12+t13;
t15 = 1.0./t14;
fu_resub = [0.0;0.0;K_tau.*t10.*t15.*1.2e+1;(K_tau.*r.*t10.*t15.*cos(x2).*1.8e+1)./Lp];
end
