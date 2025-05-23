function fx_resub = fx_jacobian(in1,u,Jrod,Jenc,Jtau,mp,r,Jtheta_corr,Lp,g,K_tau,R_tau,Ctheta,Kwire,Cenc)
parameters;
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
x4 = in1(4,:);
t2 = cos(x2);
t3 = sin(x2);
t4 = Lp.^2;
t5 = Lp.^3;
t6 = mp.^2;
t7 = r.^2;
t8 = x2.*2.0;
t9 = x2.*3.0;
t10 = x3.^2;
t11 = x4.^2;
t17 = Jenc.*1.2e+1;
t18 = Jenc.*2.4e+1;
t19 = Jrod.*1.2e+1;
t20 = Jrod.*2.4e+1;
t21 = Jtau.*1.2e+1;
t22 = Jtau.*2.4e+1;
t23 = Jtheta_corr.*1.2e+1;
t24 = Jtheta_corr.*2.4e+1;
t25 = 1.0./Lp;
t27 = 1.0./R_tau;
t12 = cos(t8);
t13 = cos(t9);
t14 = sin(t8);
t15 = t3.^2;
t16 = t3.^3;
t26 = 1.0./t4;
t28 = t4.*4.0;
t29 = t7.*9.0;
t31 = mp.*t7.*3.0;
t32 = mp.*t7.*1.5e+1;
t37 = mp.*r.*t3.*t4.*x4.*1.2e+1;
t30 = mp.*t28;
t33 = -t16;
t38 = mp.*t4.*t12.*-4.0;
t40 = mp.*t15.*t29;
t41 = t28+t29;
t42 = mp.*t7.*t12.*-9.0;
t35 = t15.*t30;
t36 = t3+t33;
t46 = t18+t20+t22+t24+t30+t32+t38+t42;
t43 = t17+t19+t21+t23+t31+t35+t40;
t47 = 1.0./t46;
t44 = 1.0./t43;
t45 = t44.^2;
et1 = Cenc.*Jenc.*R_tau.*x4.*-3.6e+1-Cenc.*Jrod.*R_tau.*x4.*3.6e+1-Cenc.*Jtau.*R_tau.*x4.*3.6e+1-Cenc.*Jtheta_corr.*R_tau.*x4.*3.6e+1-Cenc.*R_tau.*mp.*t7.*x4.*3.6e+1+R_tau.*g.*t5.*t6.*t16.*6.0+R_tau.*t2.*t6.*t10.*t16.*1.0./t26.^2.*4.0+Jenc.*Lp.*R_tau.*g.*mp.*t3.*1.8e+1+Jrod.*Lp.*R_tau.*g.*mp.*t3.*1.8e+1+Jtau.*Lp.*R_tau.*g.*mp.*t3.*1.8e+1+Jtheta_corr.*Lp.*R_tau.*g.*mp.*t3.*1.8e+1+K_tau.*Lp.*mp.*r.*t2.*u.*1.8e+1+Lp.*R_tau.*g.*t3.*t6.*t7.*1.8e+1-Cenc.*R_tau.*mp.*t4.*t15.*x4.*1.2e+1+Jenc.*R_tau.*mp.*t4.*t10.*t14.*6.0+Jrod.*R_tau.*mp.*t4.*t10.*t14.*6.0+Jtau.*R_tau.*mp.*t4.*t10.*t14.*6.0+Jtheta_corr.*R_tau.*mp.*t4.*t10.*t14.*6.0+R_tau.*t4.*t6.*t7.*t10.*t14.*6.0-R_tau.*t4.*t6.*t7.*t11.*t14.*(9.0./2.0)-Ctheta.*Lp.*R_tau.*mp.*r.*t2.*x3.*1.8e+1-Kwire.*Lp.*R_tau.*mp.*r.*t2.*x1.*1.8e+1;
et2 = R_tau.*r.*t3.*t5.*t6.*x3.*x4.*(t15-1.0).*1.2e+1;
mt1 = [0.0,0.0,Kwire.*t44.*-1.2e+1,Kwire.*r.*t2.*t25.*t44.*-1.8e+1,0.0,0.0,t25.*t47.*(Cenc.*r.*t3.*x4.*3.6e+1+Lp.*g.*mp.*r.*t12.*1.8e+1+mp.*r.*t2.*t4.*t10.*3.0-mp.*r.*t2.*t4.*t11.*1.2e+1+mp.*r.*t4.*t10.*t13.*9.0-mp.*t5.*t12.*x3.*x4.*1.6e+1)+mp.*t14.*t25.*t27.*t41.*t45.*(K_tau.*Lp.*u.*-1.2e+1+Ctheta.*Lp.*R_tau.*x3.*1.2e+1+Kwire.*Lp.*R_tau.*x1.*1.2e+1+Cenc.*R_tau.*r.*t2.*x4.*1.8e+1-Lp.*R_tau.*g.*mp.*r.*t14.*(9.0./2.0)+R_tau.*mp.*r.*t3.*t4.*t11.*6.0-R_tau.*mp.*r.*t4.*t10.*t36.*6.0+R_tau.*mp.*t5.*t14.*x3.*x4.*4.0)];
mt2 = [t25.*t27.*t47.*(K_tau.*r.*t3.*u.*-3.6e+1+Jenc.*R_tau.*g.*t2.*3.6e+1+Jrod.*R_tau.*g.*t2.*3.6e+1+Jtau.*R_tau.*g.*t2.*3.6e+1+Jtheta_corr.*R_tau.*g.*t2.*3.6e+1-Cenc.*Lp.*R_tau.*t14.*x4.*2.4e+1-R_tau.*mp.*t5.*t10.*cos(x2.*4.0).*4.0+Ctheta.*R_tau.*r.*t3.*x3.*3.6e+1+Lp.*R_tau.*t10.*t12.*t18+Lp.*R_tau.*t10.*t12.*t20+Lp.*R_tau.*t10.*t12.*t22+Lp.*R_tau.*t10.*t12.*t24+Kwire.*R_tau.*r.*t3.*x1.*3.6e+1+R_tau.*g.*mp.*t2.*t4.*9.0+R_tau.*g.*mp.*t2.*t7.*3.6e+1-R_tau.*g.*mp.*t4.*t13.*9.0+R_tau.*mp.*t5.*t10.*t12.*4.0+Lp.*R_tau.*mp.*t7.*t10.*t12.*2.4e+1-Lp.*R_tau.*mp.*t7.*t11.*t12.*1.8e+1-R_tau.*mp.*r.*t2.*t4.*x3.*x4.*6.0-R_tau.*mp.*r.*t4.*t13.*x3.*x4.*1.8e+1)-t14.*t26.*t27.*t41.*t45.*(et1+et2),1.0,0.0];
mt3 = [-t44.*(Ctheta.*1.2e+1+t14.*t30.*x4-Lp.*mp.*r.*t36.*x3.*1.2e+1),t25.*t44.*(-t37-Ctheta.*r.*t2.*1.8e+1+Lp.*t14.*t17.*x3+Lp.*t14.*t19.*x3+Lp.*t14.*t21.*x3+Lp.*t14.*t23.*x3+Lp.*mp.*t7.*t14.*x3.*1.2e+1+mp.*r.*t4.*t16.*x4.*1.2e+1+mp.*t2.*t5.*t16.*x3.*8.0),0.0,1.0,-t25.*t44.*(t37+Cenc.*r.*t2.*1.8e+1+mp.*t2.*t3.*t5.*x3.*8.0),(t26.*t44.*(Cenc.*Jenc.*6.0+Cenc.*Jrod.*6.0+Cenc.*Jtau.*6.0+Cenc.*Jtheta_corr.*6.0+Cenc.*mp.*t7.*6.0+Cenc.*mp.*t4.*t15.*2.0+r.*t3.*t5.*t6.*x3.*2.0-r.*t5.*t6.*t16.*x3.*2.0+t4.*t6.*t7.*t14.*x4.*(3.0./2.0)).*-6.0)./mp];
fx_resub = reshape([mt1,mt2,mt3],4,4);
end
