%% Optimal Control for Shuttle Dynamics - Direct Methods
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
u_traj = u_guess * (10^3) * ones(1,N);
u_norm = 300e3;                 % user defined maximum thrust for normalization

% Weightening matrices
R = 0.1/u_norm^2;               % Control input weighting matrix
p_norm = [1/x_f(1); 1/x_f(2); 1/x_i(3); 1/(2*pi)].^2;
p_weights = [10000; 100000000; 1; 10000];
P = diag(p_norm.*p_weights);    % Final state error weighting matrix
Q = zeros(4);                   % State weighting matrix (optional)

%% Dynamical system
% State dynamics
dx = @(x,u)[x(2) * sin(x(4));  
              u/x(3) - Cd*x(2)^2/x(3) - g*sin(x(4));  
              -alpha * u;       
              Cl * x(2) / x(3) - g * cos(x(4)) / x(2)];

% Jacobian of dynamics wrt state space x 
% [df1/dX(1), df1/dX(2), ...; df2/dX(1), df2/dX(2), ...; ...]
dfdx = @(x,u) [ 0,                    sin(x(4)),                               0,     x(2)*cos(x(4));
              0,              -2*Cd*x(2)/x(3),  Cd*x(2)^2/(x(3)^2) - u/(x(3)^2),       -g*cos(x(4));
              0,                            0,                               0,                  0;
              0, Cl/x(3) + g*cos(x(4))/x(2)^2,                 -Cl*x(2)/x(3)^2,   g*sin(x(4))/x(2)];

% Jacobian of dynamics wrt contro input u
dfdu = @(x,u) [ 0;
              1/x(3);
              -alpha;
              0    ];

% Running Cost function and its gradients
L = @(x,u) 0.5*x'*Q*x + 0.5*R*u.^2 ;    % Running cost
dLdx =@(x,u) Q*x;
dLdu =@(x,u) R*u;

% Final Cost function and its gradients
phi = @(x) 0.5*(x - x_f)'*P*(x - x_f) ; % Final cost
dphidx =@(x) P*(x - x_f);

%% Options section for optimum problem setup
options=optimoptions('fmincon','SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true, ...
    'FiniteDifferenceType','central', ...
    'MaxFunctionEvaluations',1e10, ...
    'FiniteDifferenceStepSize', 1e-10,...
    'OptimalityTolerance', 1e-9, ...
    'StepTolerance', 1e-9,...
    'UseParallel',true, ...
    'MaxIterations',1e3,...
    'PlotFcn', 'optimplotfval','Display','iter');

%% Setup minimization problem

nx = 4;                  % Number of states
nu = 1;                  % number of inputs

%__________Define initial conditions for the minimization problem__________

% Initialize state trajectory with initial condition
x_traj = zeros(nx, N);
% Set initial state
x_traj(:,1) = x_i; 
% Integrate the state trajectory with forward eulero
for ii = 1:N-1
    x_traj(:,ii+1) = x_traj(:,ii) + dt * dx(x_traj(:,ii), u_traj(:,ii));
end

% Combine x and u into initial guess vector z0
z0 = zeros(N*(nx + nu) + nx, 1);
for ii = 0:N-1
    % Set state guess
    z0((1 + ii*(nu + nx)):(nx + ii*(nu + nx))) = x_traj(:,ii+1);
    % Set control guess
    z0((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx))) = u_traj(:,ii+1);
end

% Add the final state to z0
z0(end-nx+1:end) = x_traj(:,end);

%___ store the necessary stuff in a structure ____________________________%

param.N = N;
param.nu = nu;
param.nx = nx;
param.dx = dx;
param.x_i = x_i;
param.dfdx = dfdx;
param.dfdu = dfdu;
param.L = L;
param.dLdx = dLdx;
param.dLdu = dLdu;
param.phi = phi;
param.dphidx = dphidx;
param.dt = dt;


%___ Define objective function and constraints ___________________________%

% define objective function
objfun = @(z) cost_and_grad(z,param);

% define non linear constraints
nonlincon = @(z) con_and_grad(z,param);

% linear inequality constraints (none in this case)
A = [];
B = [];

% linear equality constraints (none in this case)
Aeq = [];
Beq = []; 

% Lower and Upper bounds for my vector z of unkonows of dimension N*(nx+nu) + nx
h_min = -inf;          % Minimum altitude
h_max = 2e4;           % Maximum altitude (final altitude goal + overshoot)
v_min = -inf;          % Minimum velocity
v_max = 10000;         % Maximum velocity (assumption based on limits)
m_min = 10000;         % Minimum mass (cannot go under initial mass - capacità carburante)
m_max = inf;           % Initial mass limit (upper bound)
gamma_min = -inf;      % Minimum flight angle
gamma_max = inf;       % Maximum flight angle
u_min = 0;             % Minimum thrust
u_max = 600e3;         % Maximum thrust

% Initialize lower and upper bounds for z
LB = zeros(size(z0));  % Lower bound vector
UB = zeros(size(z0));  % Upper bound vector

% Set bounds for each time step in z
for i = 0:N-1
    % State bounds (altitude, velocity, mass, flight angle)
    LB((1 + i * (nx + nu)):(nx + i * (nx + nu))) = [h_min; v_min; m_min; gamma_min];
    UB((1 + i * (nx + nu)):(nx + i * (nx + nu))) = [h_max; v_max; m_max; gamma_max];
    
    % Control bounds (thrust)
    LB((1 + nx + i * (nx + nu)):(nx + nu + i * (nx + nu))) = u_min;
    UB((1 + nx + i * (nx + nu)):(nx + nu + i * (nx + nu))) = u_max;
end

% Final state bounds
LB(end-nx+1:end) = [h_min; v_min; m_min; gamma_min];
UB(end-nx+1:end) = [h_max; v_max; m_max; gamma_max];


%% CheckGradients
check_options=optimoptions("fmincon",FiniteDifferenceType="central");

[valid_cost,err_gradcost]=checkGradients(@(z) cost_and_grad(z,param), z0,check_options,'Display','on');
[valid_con,err_gradcon]=checkGradients(@(z) con_and_grad(z,param), z0,check_options,'IsConstraint',true,'Display','on');

disp('CheckGradients:')
disp(['- Functional Cost: ',num2str(valid_cost)]);
disp(['- Inequality Constraints: ',num2str(valid_con(1))]);
disp(['- Equality Constraints: ',num2str(valid_con(2))]);



%% Simulation run

bounds = input('Scegli se usare i bounds o no (Y=1 N=0) : ');
if bounds
    tic
    [z,fval] = fmincon(objfun,z0,A,B,Aeq,Beq,LB,UB,nonlincon,options);
else
    tic
    [z,fval] = fmincon(objfun,z0,A,B,Aeq,Beq,[],[],nonlincon,options);
end

tempo = toc;
disp(['Optimization Time: ', num2str(toc), ' sec.']);
disp(['Final Functional Cost: J=', num2str(fval)]);

%% Final plots
figure
subplot(2,2,1)
plot(Tu,z(1:(nx+nu):N*(nx+nu)-nx),'k','LineWidth',2); hold on; grid on;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Altitude [m]', 'FontSize', 14, 'Interpreter', 'latex');
yline(x_f(1), 'k--', 'LineWidth', 1)
legend('Altitude', '$h_r$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');

subplot(2,2,2)
plot(Tu,z(2:(nx+nu):N*(nx+nu)-nx+2),'b','LineWidth',2); hold on; grid on;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Velocity [m/s]', 'FontSize', 14, 'Interpreter', 'latex');
yline(x_f(2), 'k--', 'LineWidth', 1)
legend('Velocity', '$v_r$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');

subplot(2,2,3)
plot(Tu,z(3:(nx+nu):N*(nx+nu)-nx+3),'r','LineWidth',2); hold on; grid on;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Mass [Kg]', 'FontSize', 14, 'Interpreter', 'latex');
yline(x_f(3), 'k--', 'LineWidth', 1)
legend('Mass','$m_r$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');

subplot(2,2,4)
plot(Tu,z(4:(nx+nu):N*(nx+nu)-nx+4),'m','LineWidth',2); hold on; grid on;
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Flight Angle [rad]', 'FontSize', 14, 'Interpreter', 'latex');
yline(x_f(4), 'k--', 'LineWidth', 1)
legend('Flight Angle', '$\gamma_r$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best');

figure
yline(u_guess*10^3, 'k--', 'LineWidth', 1)
hold on
plot(Tu(1:end-1),z(5:(nx+nu):N*(nx+nu)-nx),'LineWidth',2);
xlabel('Time [s]', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('Thrust [N]', 'FontSize', 14, 'Interpreter', 'latex');
title('Control Input (Thrust)', 'FontSize', 14, 'Interpreter', 'latex');
legend('u first guess', '$u_{optimal}$', 'FontSize', 14, 'Interpreter', 'latex', 'Location', 'best')
axis tight

%% COST/OBJECTIVE function {J} AND GRADIENT [grad(J)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cost,grad] = cost_and_grad(z,param)

N = param.N;
nx = param.nx;
nu = param.nu;
dt = param.dt;
L = param.L;
phi = param.phi;

% extract states and control from z
x = zeros(nx,N+1); u = zeros(nu,N);
for ii = 0:N
    x(:,ii+1) = z((1 + ii*(nu + nx)):(nx + ii*(nu + nx)));
end    
for ii = 0:N-1
    u(:,ii+1) = z((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx)));    
end
    
% Cost Function 
cost = 0;
for ii = 1:N
    cost = cost + dt*L(x(:,ii),u(:,ii));
end
cost = cost + phi(x(:,end));

if nargout > 1
% Gradient of the objective function
dLdx = param.dLdx;
dLdu = param.dLdu;
dphidx = param.dphidx;

grad = zeros(size(z));
for ii = 0:N-1
    grad((1 + ii*(nu + nx)):(nx + ii*(nu + nx)),1) = dt*dLdx(x(:,ii + 1),u(:,ii + 1));    
    grad((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx)),1) = dt*dLdu(x(:,ii + 1),u(:,ii + 1));    
end
grad(end - nx + 1:end,1) = dphidx(x(:,end));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CONSTRAINTS AND GRADIENTS function
function [c,con,g,grad] = con_and_grad(z,param)

% extract parameters

N = param.N;
nx = param.nx;
nu = param.nu;
dx = param.dx;
dt = param.dt;
x_i = param.x_i;

% extract states and control from z

x = zeros(nx,N+1); u = zeros(nu,N);
for ii = 0:N
    x(:,ii+1) = z((1 + ii*(nu + nx)):(nx + ii*(nu + nx)));
end    
for ii = 0:N-1
    u(:,ii+1) = z((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx)));    
end

% constraint function

c = []; % here is room for inequality constraint
con = zeros(N*nx,1);
con(1:nx) = x_i - x(:,1) ; % initial condition constraint

for ii = 1:N
    con((1:nx) + ii*nx) = x(:,ii) + dt*dx(x(:,ii),u(:,ii)) - x(:,ii+1);    
end

% gradient of the constraint

dfdx = param.dfdx;
dfdu = param.dfdu;

if nargout > 2

    g = []; % here is room for inequality constraint
    grad = zeros(nx,N*(nu+nx) + nx); % gradient of a vector
    grad(1:nx,1:nx) = - eye(nx); 
    for ii = 1:N
        grad((1 + nx +(ii - 1)*nx):((nx + ii*nx)),(1 +(ii - 1)*(nx+nu)):((ii - 1)*(nx+nu) + nx)) = ...
             eye(nx) + dt*dfdx(x(:,ii),u(:,ii));
        grad((1 + nx +(ii - 1)*nx):((nx + ii*nx)),(1 +(ii)*(nx+nu)):((ii)*(nx+nu) + nx)) = ...
             - eye(nx);
        grad((1 + nx +(ii - 1)*nx):((nx + ii*nx)),(1 +(ii - 1)*(nx+nu) + nx):((ii - 1)*(nx+nu) + nx + nu)) = ...
             + dt*dfdu(x(:,ii),u(:,ii));
    end
    grad = grad.';
    
end

end

