%% Optimal Control for Qube Double Pendulum Swing Up - Direct Methods
% LAB - Mechatronics A AA.24/25
% Lorenzo Bauce ** lorenzo.bauce@mail.polimi.it **
% All rights reserved

close all;
clear;
clc;

% Initial state : theta0, phi0, dtheta0, dphi0
x_i = [0; -pi; 0; 0];

% Final state : theta_f, phi_f, dtheta_f, dphi_f
x_f = [0; 0; 0; 0];

% Time settings
t0 = 0;
tf = 1; 
ts = 0.002;     % Sampling time of the QUBE Servo 2 micro-controller
dt = input('Time step for OCP: ');
Tu = t0:dt:tf;
Ns = tf/ts +1;
N = tf/dt + 1;

% Control input bounds
u_max = 10; % volt
u_min = -10; % volt

% Initial guess for control input history
u_traj = ones(1,N);

% Weighting matrices
p_weight = [0.01/(2*pi)^2, 1000/(2*pi)^2, 0.01/(20)^2, 0.1/(20)^2]; % normalized weighting vector
P = diag(p_weight);  % weight on final state error
Q = zeros(4);        % weight on state
R = 1/u_max^2;       % weight on control input



%% Dynamics and Jacobians
% import parameters (to be updated)
parameters;

% State dynamics
dx= @(x, u) Eq_pend_inv(x,u);

% Jacobian of f wrt x
dfdx = @(x,u) fx_jacobian(x,u);

% Jacobian of f wrt u
dfdu = @(x,u) fu_jacobian(x,u);

% Running Cost function and its gradients
L = @(x,u) 0.5*x'*Q*x + 0.5*R*u.^2 ; % Running cost
dLdx =@(x,u) Q*x;
dLdu =@(x,u) R*u;

% Final Cost function and its gradients
phi = @(x) 0.5*(x - x_f)'*P*(x - x_f) ; % Final cost
dphidx =@(x) P*(x - x_f);

%% Options section for optimum problem setup
options=optimoptions( ...
    'fmincon','SpecifyObjectiveGradient',true, ...
    'Algorithm','interior-point',...
    'ConstraintTolerance', 1e-8, ...
    'SpecifyConstraintGradient', true, ...
    'FiniteDifferenceType','central', ...
    'MaxFunctionEvaluations',1e10, ...
    'FiniteDifferenceStepSize', 1e-4,...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-7,...
    'UseParallel',true, ...
    'MaxIterations',1e3,...
    'PlotFcn', 'optimplotfval','Display','iter');


%% Setup minimization problem
nx = 4;    % Number of states
nu = 1;    % Number of inputs

% Ensure that the number of intervals is N. Here, we have N+1 state points.
x_traj = zeros(nx, N+1);
x_traj(:,1) = x_i;  % Set the initial condition

% Generate the state trajectory using Simpson's rule integration
% (Assumes that u_traj is provided for each interval, i.e., for ii = 1,...,N)
for ii = 1:N
    % f1: dynamics at the beginning of the interval
    f1 = dx(x_traj(:,ii), u_traj(:,ii));
    
    % Estimate the endpoint using an Euler step
    x_temp = x_traj(:,ii) + dt * f1;
    
    % f3: dynamics at the estimated endpoint
    f3 = dx(x_temp, u_traj(:,ii));
    
    % Estimate the midpoint state (using Euler half-step)
    x_mid = x_traj(:,ii) + (dt/2) * f1;
    
    % f2: dynamics at the midpoint
    f2 = dx(x_mid, u_traj(:,ii));
    
    % Simpson's rule update for the state
    x_traj(:,ii+1) = x_traj(:,ii) + (dt/6) * (f1 + 4*f2 + f3);
end

% Combine state and control trajectories into the decision variable vector z0.
% The structure is assumed to be:
%    z0 = [ x(:,1); u(:,1); x(:,2); u(:,2); ...; x(:,N); u(:,N); x(:,N+1) ]
z0 = zeros(N*(nx + nu) + nx, 1);
for ii = 0:N-1
    % Insert state guess for time step ii+1
    z0((1 + ii*(nu + nx)):(nx + ii*(nu + nx))) = x_traj(:,ii+1);
    
    % Insert control guess for time step ii+1
    z0((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx))) = u_traj(:,ii+1);
end


% Add the final state to z0
z0(end-nx+1:end) = x_traj(:,end);


%---- store the necessary stuff in a structure ----%

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

%---- define objective function and constraints ---%

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
theta_min = -pi/2;
theta_max = +pi/2;
phi_min = -inf;
phi_max = +inf;
dtheta_min = -inf;
dtheta_max = +inf;
dphi_min = -inf;
dphi_max = +inf; 

% Initialize lower and upper bounds for z
LB = zeros(size(z0));  % Lower bound vector
UB = zeros(size(z0));  % Upper bound vector

% Set bounds for each time step in z
for i = 0:N-1

    % State bounds
    LB((1 + i * (nx + nu)):(nx + i * (nx + nu))) = [theta_min; phi_min; dtheta_min; dphi_min];
    UB((1 + i * (nx + nu)):(nx + i * (nx + nu))) = [theta_max; phi_max; dtheta_max; dphi_max];
    
    % Control bounds
    LB((1 + nx + i * (nx + nu)):(nx + nu + i * (nx + nu))) = u_min;
    UB((1 + nx + i * (nx + nu)):(nx + nu + i * (nx + nu))) = u_max;
end

% Final state bounds
LB(end-nx+1:end) = [theta_min; phi_min; dtheta_min; dphi_min];
UB(end-nx+1:end) = [theta_max; phi_max; dtheta_max; dphi_max];

%% CheckGradients
check_options=optimoptions("fmincon",FiniteDifferenceType="central");

[valid_cost,err_gradcost]=checkGradients(@(z) cost_and_grad(z,param), z0,check_options,'Display','on');
[valid_con,err_gradcon]=checkGradients(@(z) con_and_grad(z,param), z0,check_options,'IsConstraint',true,'Display','on');

disp('CheckGradients:')
disp(['- Functional Cost: ',num2str(valid_cost)]);
disp(['- Inequality Constraints: ',num2str(valid_con(1))]);
disp(['- Equality Constraints: ',num2str(valid_con(2))]);


%% Simulation run
bounds = input('Do you want to use state bounds? [(Y=1 N=0)] : ');
if bounds
    tic
    [z,fval] = fmincon(objfun,z0,A,B,Aeq,Beq,LB,UB,nonlincon,options);
else
    tic
    [z,fval] = fmincon(objfun,z0,A,B,Aeq,Beq,[],[],nonlincon,options);
end

toc;
%disp(['Optimization Time is ', num2str(toc), ' seconds']);
disp(['Final Functional Cost is ', num2str(fval)]);

%% Final Plots
figure("Name",'QUBE Dynamics with Optimal Control');
subplot(2,1,1)
hold on
grid on
for kk = 1:2
    if kk == 1
        plot(Tu(1:end),z(kk:(nx+nu):N*(nx+nu)-nx),'LineWidth',1);
    else
        plot(Tu(1:end),z(kk:(nx+nu):N*(nx+nu)-nx+kk),'LineWidth',1);

    end
end
yline(x_f(2),'k--','LineWidth',1.5)
yline(20*pi/180, 'r-.')
yline(-20*pi/180, 'r.-')
legend('\theta', '\phi', '\phi ref', 'control switch', 'Location','best')
axis tight

subplot(2,1,2)
hold on
grid on
for kk = 3:4
    if kk == 1
        plot(Tu(1:end),z(kk:(nx+nu):N*(nx+nu)-nx),'LineWidth',1);
    else
        plot(Tu(1:end),z(kk:(nx+nu):N*(nx+nu)-nx+kk),'LineWidth',1);

    end
end
legend('\delta\theta', '\delta \phi', 'Location','best')
axis tight

figure
plot(Tu(1:end),z(5:(nx+nu):N*(nx+nu)-nx+5),'LineWidth',1);
legend('Optimal Control u(t)')
grid on
ylim([-10 10])

%% Dynamics solved by accurate integrator [ode45]
% Extract optimal control trajectory
U_opt = z(5:(nx + nu):N * (nx + nu) - nx + 5);
tempo = 0:0.002:tf;
U_opt_q = interp1(Tu,U_opt,tempo);

% Initial state
x0 = x_i;

% Solve dynamics using ode45
[T,X]=ode45(@(t,x) stateEq_pend_inv(t,x,U_opt_q,tempo),tempo,x0);

% Plot results
figure("Name",'QUBE Dynamics OC integrated with ode45');
subplot(2,1,1)
plot(T, X(:,1),'m',T, X(:,2),'c', 'LineWidth',1);
xlabel('Time (s)');
ylabel('Angle (rad)');
legend('\theta', '\phi');
grid on;

subplot(2,1,2)
plot(T, X(:,3),'m',T, X(:,4),'c','LineWidth',1);
xlabel('Time (s)');
ylabel('Angular velocity (rad/s)');
legend('d\theta', 'd\phi');
title('Pendulum Dynamics with Optimal Control');
grid on;

%% COST function {J} AND GRADIENT {grad(J)}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cost,grad] = cost_and_grad(z,param)

% Parameters from param
N   = param.N;    % number of control points (assumed odd so that N-1 is even)
nx  = param.nx;
nu  = param.nu;
dt  = param.dt;
L   = param.L;    % running cost function: L(x,u)
phi = param.phi;  % terminal cost function: phi(x)

% Extract states and controls from z
% x has N+1 points, u has N points
x = zeros(nx, N+1); 
u = zeros(nu, N);
for ii = 0:N
    % x(:,1) corresponds to time 0, x(:,N+1) is the final state
    x(:,ii+1) = z((1 + ii*(nu + nx)):(nx + ii*(nu + nx)));
end    
for ii = 0:N-1
    % u(:,1) corresponds to time 0, u(:,N) is the last control input
    u(:,ii+1) = z((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx)));    
end
    
% --- Cost Function using Simpson's rule ---
% We approximate the integral of L(x,u) over the horizon corresponding to
% the N control points. (Assumes N is odd so that the number of subintervals, N-1, is even.)
cost_running = 0;
for ii = 1:N
    % Simpson coefficient:
    % weight = 1 at the endpoints, 4 for even indices (interior) and 2 for odd interior indices.
    if ii == 1 || ii == N
        weight = 1;
    elseif mod(ii,2) == 0
        weight = 4;
    else
        weight = 2;
    end
    cost_running = cost_running + weight * L(x(:,ii), u(:,ii));
end
% Multiply by dt/3 per Simpson's rule and add terminal cost.
cost = (dt/3) * cost_running + phi(x(:,end));

% --- Gradient of the objective function ---
if nargout > 1
    dLdx   = param.dLdx;
    dLdu   = param.dLdu;
    dphidx = param.dphidx;

    grad = zeros(size(z));
    for ii = 0:N-1
        index = ii + 1;  % corresponds to the sample index in 1,...,N
        if index == 1 || index == N
            weight = 1;
        elseif mod(index,2) == 0
            weight = 4;
        else
            weight = 2;
        end
        % Gradient with respect to the state at time index "index"
        grad((1 + ii*(nu + nx)):(nx + ii*(nu + nx)),1) = (dt/3) * weight * dLdx(x(:,index), u(:,index));    
        % Gradient with respect to the control at time index "index"
        grad((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx)),1) = (dt/3) * weight * dLdu(x(:,index), u(:,index));    
    end
    % Add terminal cost gradient (with respect to the final state x(:,N+1))
    grad(end - nx + 1:end,1) = dphidx(x(:,end));
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dynamic constraints f(x(t),u(t),t) and gradients {df/dx, df/du}
function [c,con,g,grad] = con_and_grad(z,param)

% Extract parameters
N   = param.N;    % number of intervals
nx  = param.nx;   % state dimension
nu  = param.nu;   % control dimension
dx  = param.dx;   % dynamics function: dx = dx(x,u)
dt  = param.dt;   % time step
x_i = param.x_i;  % given initial state

% Extract states and controls from z
% States: x(:,1) ... x(:,N+1); Controls: u(:,1) ... u(:,N)
x = zeros(nx, N+1);
u = zeros(nu, N);
for ii = 0:N
    x(:,ii+1) = z((1 + ii*(nu + nx)):(nx + ii*(nu + nx)));
end    
for ii = 0:N-1
    u(:,ii+1) = z((1 + nx + ii*(nu + nx)):(nx + nu + ii*(nu + nx)));    
end

% Constraint function
%
% We now have two sets of constraints:
% 1. Initial condition constraint: x(:,1) == x_i.
% 2. Dynamics constraints using Simpson's rule.
%
% We assemble these into one vector, con. (c is left empty here for inequality constraints.)
c = [];  % room for inequality constraints (none here)
% Total constraints: initial condition (nx) + N intervals (N*nx)
con = zeros((N+1)*nx, 1);

% (a) Initial condition constraint
con(1:nx) = x_i - x(:,1);

% (b) Dynamics constraints for each interval using Simpson's rule
for ii = 1:N
    % For the interval from x(:,ii) to x(:,ii+1) with control u(:,ii):
    % Approximate the integral using Simpson's rule.
    % Use the midpoint approximation: x_mid ~ (x(:,ii) + x(:,ii+1))/2.
    x_mid = (x(:,ii) + x(:,ii+1)) / 2;
    
    % Simpson's integration of the dynamics:
    % x(:,ii) + (dt/6)*( dx(x(:,ii), u(:,ii)) + 4*dx(x_mid, u(:,ii)) + dx(x(:,ii+1), u(:,ii)) )
    % should equal x(:,ii+1). Thus, we set:
    con((1:nx) + ii*nx) = x(:,ii) + (dt/6)*( dx(x(:,ii), u(:,ii)) + 4*dx(x_mid, u(:,ii)) + dx(x(:,ii+1), u(:,ii)) ) - x(:,ii+1);
end

% Gradient of the constraints
%
% We assume that param contains the functions:
%   dfdx(x,u) = partial derivative of dx with respect to x,
%   dfdu(x,u) = partial derivative of dx with respect to u.
%
% The structure of the decision vector z is as before:
%   [ x(:,1); u(:,1); x(:,2); u(:,2); ... ; x(:,N); u(:,N); x(:,N+1) ].
%
% The gradient "grad" will be built so that each dynamic constraint
% row depends on x(:,ii), u(:,ii), and x(:,ii+1) as follows.
dfdx_fun = param.dfdx;
dfdu_fun = param.dfdu;

if nargout > 2
    g = [];  % no inequality constraint gradients here
    
    % Total number of decision variables:
    % states: (N+1)*nx, controls: N*nu
    nvar = (N+1)*nx + N*nu;
    % Total number of constraints: (N+1)*nx.
    % We form a matrix "grad" of size (constraints) x (nvar) and then transpose it.
    grad = zeros((N+1)*nx, nvar);
    
    % --- Initial condition constraint: con(1:nx) = x_i - x(:,1)
    % Its derivative w.r.t. x(:,1) is -I.
    grad(1:nx, 1:nx) = -eye(nx);
    
    % --- Dynamic constraints for each interval (ii = 1,...,N)
    for ii = 1:N
        % The constraint for interval ii is associated with rows:
        row_idx = (1:nx) + ii*nx;
        
        % Identify the decision variable indices:
        % For x(:,ii): it is stored in z starting at index:
        col_idx_xk = (1 + (ii-1)*(nx+nu)) : ((ii-1)*(nx+nu) + nx);
        % For u(:,ii):
        col_idx_uk = (1 + nx + (ii-1)*(nx+nu)) : ((ii-1)*(nx+nu) + nx + nu);
        % For x(:,ii+1):
        col_idx_xkp1 = (1 + ii*(nx+nu)) : (ii*(nx+nu) + nx);
        
        % Extract the states and control for this interval:
        xk   = x(:,ii);
        xkp1 = x(:,ii+1);
        uk   = u(:,ii);
        % The midpoint is approximated as:
        x_mid = (xk + xkp1)/2;
        
        % Evaluate the Jacobians of the dynamics at the three sample points:
        A   = dfdx_fun(xk,   uk);      % at xk
        B   = dfdx_fun(x_mid, uk);     % at the midpoint (chain rule gives a factor 1/2 later)
        C   = dfdx_fun(xkp1, uk);      % at xkp1
        
        A_u = dfdu_fun(xk,   uk);       % control derivative at xk
        B_u = dfdu_fun(x_mid, uk);      % control derivative at midpoint
        C_u = dfdu_fun(xkp1, uk);       % control derivative at xkp1
        
        % The defect constraint is:
        % F = xk + (dt/6)[ f(xk,uk) + 4 f(x_mid,uk) + f(xkp1,uk) ] - xkp1 = 0.
        %
        % Its partial derivatives are computed as follows:
        %
        % With respect to xk:
        %   dF/dxk = I + (dt/6)*[ dfdx(xk,uk) + 4*(1/2)*dfdx(x_mid,uk) ]
        %           = I + (dt/6)*[ A + 2*B ].
        dFdxk = eye(nx) + (dt/6)*( A + 2*B );
        
        % With respect to xkp1:
        %   dF/dxkp1 = (dt/6)*[ 4*(1/2)*dfdx(x_mid,uk) + dfdx(xkp1,uk) ] - I
        %             = (dt/6)*( 2*B + C ) - eye(nx).
        dFdxkp1 = (dt/6)*( 2*B + C ) - eye(nx);
        
        % With respect to u (which appears in all three evaluations):
        %   dF/du = (dt/6)*[ dfdu(xk,uk) + 4*dfdu(x_mid,uk) + dfdu(xkp1,uk) ]
        dFdu = (dt/6)*( A_u + 4*B_u + C_u );
        
        % Insert these blocks into the gradient matrix.
        grad(row_idx, col_idx_xk)   = dFdxk;
        grad(row_idx, col_idx_uk)   = dFdu;
        grad(row_idx, col_idx_xkp1) = dFdxkp1;
    end
    
    grad = grad.';
end

end


