%% Optimal Control for Shuttle Dynamics - Indirect Methods
% Mechatronic systems A AA.24/25
% Lorenzo Bauce ** lorenzo.bauce@mail.polimi.it **
% All rights reserved

close all; clear; clc;

%% Data
g = 9.81;          % Gravitational acceleration (m/s^2)
Cd = 0.5;          % Drag coefficient
Cl = 0.2;          % Lift coefficient
alpha = 0.0002;    % Fuel consumption


x_i = [0; 1; 25400; pi/2];      % Initial state
x_f = [8000; 300; 21400; pi];   % Final state

t0 = 0;
tf = 80;                        % Finite time horizon 
dt = input('Time step? ');      % Select 0.1 to 0.5
Tu = t0:dt:tf;                  % Time span
N = tf/dt + 1;                  % Time instants

u_guess = input('Constant control Guess in [KN]? '); % 350 used for report
u = u_guess * (10^3) * ones(1,N);

R = 0.1/(max(u))^2;             % Control input weighting matrix
p_norm = [1/x_f(1); 1/x_f(2); 1/x_i(3); 1/(2*pi)].^2;
p_weights = [100000; 10000; 1; 100];
P = diag(p_norm.*p_weights);    % Final state error weighting matrix
Q = zeros(4);                   % State weighting matrix

%% Starting iterations
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6 * ones(4,1));
step_max = 1e5; step_min = 1e4;   % speed of control adjustments                        

% First order necessary condition for optimality in fixed end-time ocp is
% delta_H/delta_u = 0 on optimal trajectory
eps = 1e-5; % condition to stop iterations

iter_max = 1e4;

% Cost function and the derivative dH/du norm progression over iterations
CostF = zeros(iter_max,1);
NormH = zeros(iter_max,1);
figure
subplot(2,1,1)
h = semilogy(NaN, NaN, 'bo', 'MarkerSize',3);
xlabel('Iteration');
ylabel('Cost J');
legend('Cost Functional variation over iterations');
grid on;
hold on;
subplot(2,1,2)
z = semilogy(NaN, NaN, 'mo', 'MarkerSize',3);
xlabel('Iteration');
ylabel('Norm \deltaH/\deltau');
legend('Norm \deltaH/\deltau variation over iterations');
grid on;
hold on;


global Stop
Stop = false;
hButton = uicontrol('Style', 'pushbutton', 'String', 'Stop', ...
                    'Position', [20 20 50 20], ...
                    'Callback', 'global Stop; Stop = true;');

tic
for ii = 1:iter_max

   % 0) Select the speed of control adjustment
   step = (step_max-step_min) .* exp(3 * (ii-1)/iter_max) + step_min;
   %step = step_max;

   % 1) Start with assumed control u and integrate state forward in time
   [Tx,X] = ode45(@(t,x) stateEq(t,x,u,Tu,Cd,Cl,alpha,g),[t0 tf], x_i,options);

   % 2) Adjoint vector trajectory moving backward in time from the final time
   X_tf = X(end,:)';            % final state vector
   p_tf = P*(X_tf - x_f);       % final adjoint vector
   [Tlmb,lmb] = ode45(@(t,lmb) AdjointEq(t,lmb,u,Tu,X,Tx,Cd,Cl,g),[tf t0],p_tf,options);

   lmb1 = lmb(:,1);
   lmb1 = interp1(Tlmb,lmb1,Tx);
   lmb2 = lmb(:,2);
   lmb2 = interp1(Tlmb,lmb2,Tx);
   lmb3 = lmb(:,3);
   lmb3 = interp1(Tlmb,lmb3,Tx);
   lmb4 = lmb(:,4);
   lmb4 = interp1(Tlmb,lmb4,Tx);

   % 3) Calculate deltaH with x2(t), x3(t), lmb2(t), lmb3(t)
   dH = dHdu(lmb1,lmb2,lmb3,lmb4,Tx,X,u,Tu,R,alpha);

   % 4) Calculate the cost functional
   J(ii,1) = 0.5*(X_tf - x_f).'*P*(X_tf - x_f) + 0.5*R*(u*u')*tf/length(Tu);
    
   % 5) Condition if dH/du < epsilon, exit
   if isnan(dH)
       disp('dH non è un numero, ritenta!');
       return;
   elseif norm(dH) < eps
       break;
   else
       % 6) adjust control for next iteration
       u_old = u;
       u = IterControl(dH,Tx,u_old,Tu,step);
   end

   % 6) Plot section
    CostF(ii) = J(ii,1);
    NormH(ii) = norm(dH);
    set(h, 'XData', 1:ii, 'YData', CostF(1:ii));
    set(z, 'XData', 1:ii, 'YData', NormH(1:ii));
    if ii ==1
        drawnow;
        disp(['Cost functional at iteration ', num2str(ii-1), ' is ', num2str(J(ii,1))]);
        disp(['Norm of dH/du at iteration ', num2str(ii-1), ' is ', num2str(NormH(ii))]);
        fprintf('\n')
    elseif mod(ii,10) == 0
        drawnow;
        disp(['Cost functional at iteration ', num2str(ii), ' is ', num2str(J(ii,1))]);
        disp(['Norm of dH/du at iteration ', num2str(ii), ' is ', num2str(NormH(ii))]);
        fprintf('\n')
    end

    global Stop
    if Stop
        disp('Iteration manually stopped');
        break;
    end
end



tempo=toc;
disp(['Final iteration: ',num2str(ii)]);
disp(['Final cost: ',num2str(J(ii,1)),' [-]']);
disp(['Optimization Time: ', num2str(toc), ' sec.']);

%% Final Plots for report
figure
subplot(2,2,1)
plot(Tx,X(:,1),'k','LineWidth',2); hold on; grid on;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Altitude [m]', 'FontSize', 14, 'Interpreter', 'latex');
yline(x_f(1), 'k--', 'LineWidth', 1)
legend('Altitude', '$h_r$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');

subplot(2,2,2)
plot(Tx,X(:,2),'b','LineWidth',2); hold on; grid on;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Velocity [m/s]', 'FontSize', 14, 'Interpreter', 'latex');
yline(x_f(2), 'k--', 'LineWidth', 1)
legend('Velocity', '$v_r$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');

subplot(2,2,3)
plot(Tx,X(:,3),'r','LineWidth',2); hold on; grid on;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Mass [Kg]', 'FontSize', 14, 'Interpreter', 'latex');
yline(x_f(3), 'k--', 'LineWidth', 1)
legend('Mass','$m_r$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');

subplot(2,2,4)
plot(Tx,X(:,4),'m','LineWidth',2); hold on; grid on;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Flight Angle [rad]', 'FontSize', 14, 'Interpreter', 'latex');
yline(x_f(4), 'k--', 'LineWidth', 1)
legend('Flight Angle', '$\gamma_r$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');

figure
yline(u_guess*10^3, 'k--', 'LineWidth', 1)
hold on
plot(Tu,u,'LineWidth',2);
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Thrust [N]', 'FontSize', 14, 'Interpreter', 'latex');
title('Control Input (Thrust)', 'FontSize', 14, 'Interpreter', 'latex');
legend('u first guess', '$u_{optimal}$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best')

%% FUNCTIONS SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = stateEq(t,x,u,Tu,Cd,Cl,alpha,g)
% Interpolate the control
u = interp1(Tu,u,t);    

% Integrate the state trajectory forwards
dx = zeros(4,1);        
dx(1) = x(2) * sin(x(4));                       % Altitude rate
dx(2) = u/x(3) - Cd*x(2)^2/x(3) - g*sin(x(4));  % Velocity rate
dx(3) = -alpha * u;                             % Mass rate
dx(4) = Cl * x(2)/x(3) - g*cos(x(4))/x(2);      % Flight angle rate
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dlmb = AdjointEq(t,lmb,u,Tu,X,Tx,Cd,Cl,g)
% Interpolate the control
u = interp1(Tu,u,t);   

% Interpolate the state variables
x1 = interp1(Tx,X(:,1),t);   
x2 = interp1(Tx,X(:,2),t);
x3 = interp1(Tx,X(:,3),t);
x4 = interp1(Tx,X(:,4),t);

% Integrate the adjoint vector (Lagrange multipliers) trajectory backwards
dlmb = zeros(4,1);
dlmb(1) = 0;
dlmb(2) = -sin(x4)*lmb(1) + 2*Cd*x2/x3*lmb(2) - g*cos(x4)/(x2^2)*lmb(4) - Cl/x3*lmb(4);
dlmb(3) = (u-Cd*x2^2)/(x3^2)*lmb(2) + Cl*x2/(x3^2)*lmb(4);
dlmb(4) = -x2*cos(x4)*lmb(1) + g*cos(x4)*lmb(2) - g*sin(x4)/x2*lmb(4);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dH = dHdu(lmb1,lmb2,lmb3,lmb4,Tx,X,u,Tu,R,alpha)
% Interpolate the control
u = interp1(Tu,u,Tx);   
% Mass state variable (x3)
x3 = X(:,3);            
% Gradient of the Hamiltonian w.r.t. u
dH = R*u + (1./x3).*lmb2 - alpha*lmb3;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_new = IterControl(dH,tx,u,tu,step)
% interploate dH/du
dH = interp1(tx,dH,tu);
% adjust the control action with relation
u_new = u - step*dH;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%