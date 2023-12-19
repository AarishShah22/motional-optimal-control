clc
clear all
close all

%% Define Parameters

a = 1; % rear axle to COG (centre of gravity)
L = 2; % total length of vehicle
tf = 2; % simulation end time, simulation begins at t=0
T = 0.05; % time step
N = tf/T; % no. of time steps
init_pose_e = [1,2.2,0]; % ego vehicle initial pose (pose = [x; y; theta]) 
% reference trajectory also follows [x,y,theta]
ref_traj = [linspace(1,20,N)', 1.75*ones(N,1), zeros(N,1)]; 
pedestrian_position = [1.5,2.75];
dist_from_ped = 1;  

% Priority structure for iterative rule-relaxation algorithm
% Order of rules is: 
%   1. Maintain safe distance from pedestrian
%   2. Lane-keeping
%   3. Min. speed limit
% 1 means the rule is enforced, 0 means the rule is relaxed
priority_structure = [1,1,1;
                      1,1,0;
                      1,0,1;
                      1,0,0;
                      0,1,1;
                      0,1,0;
                      0,0,1;
                      0,0,0];

%% Rule Relaxation Algorithm

% Iterative rule-relaxation algorithm, relaxes rules according to priority
% structure and checks for feasible solution. If found, exits the loop. If
% not, continues to next scenario and checks for feasible solution.
for i = 1:length(priority_structure)
    i
    % information about which rules are enforced/relaxed
    enforced_rules = priority_structure(i,:); 
    % initial guess for fmincon
    u0 = [10*ones(N,1) ones(N,1)*0.3 ones(N,1)*10 ones(N,1)*0 ones(N,1)*1 ones(N,1)*0];
    % using SQP algorithm with fmincon since we have non-linear dynamics
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', ...
        'sqp','StepTolerance',1e-6, 'MaxFunctionEvaluations',1e6, ...
        'MaxIterations',1e3);
    [u_opt, fval, exit_flag] = fmincon( ...
                               @(u) objective_fun(u), u0, [], [], [], ...
                               [], [], [],@(u) nl_constraints(u, ...
                               init_pose_e,ref_traj, pedestrian_position, ...
                               N,T,enforced_rules, dist_from_ped), options);
    % if optimization is successful and feasible solution is found, exit
    % the for loop and return that solution. Else continue on to next
    % scenario.
    if exit_flag == 1 || exit_flag == 2
        break;
    end
end

% Construct state trajectory based on optimal control trajectory
traj_opt = zeros(N,3);
traj_opt(1,1) = init_pose_e(1);
traj_opt(1,2) = init_pose_e(2);
traj_opt(1,3) = init_pose_e(3);

for i =1:N
    traj_opt(i+1,1) = traj_opt(i,1) + u_opt(i,1)*cos(traj_opt(i,3))*T;
    traj_opt(i+1,2) = traj_opt(i,2) + u_opt(i,1)*sin(traj_opt(i,3))*T;
    traj_opt(i+1,3) = traj_opt(i,3) + u_opt(i,2)*T;
end

%% Plotting

t = linspace(0,tf,N+1);

% Plotting optimal trajectory
figure(1)
plot(traj_opt(:,1),traj_opt(:,2))
hold on
scatter(traj_opt(:,1),traj_opt(:,2))
axis square
axis equal
xlim([0, 20]); % Set x-axis limits
ylim([0, 6]); % Set y-axis limits
plot([0,20],[2,2],'Linewidth',2)
plot([0,20],[0.75,0.75],'Linewidth',2)
% draw the safe circle around the pedestrian
th = 0:pi/50:2*pi;
xunit = dist_from_ped*cos(th) + pedestrian_position(1);
yunit = dist_from_ped*sin(th) + pedestrian_position(2);
plot(xunit, yunit);
xlabel("X")
ylabel("Y")
legend("Trajectory","","Relaxable Lane Limit", "Unrelaxable Lane Limit", "Safe Circle Around Pedestrian")
title('Trajectory Generated')

% Plotting states v/s time
figure(2)
plot(t,traj_opt(:,1))
hold on
plot(t,traj_opt(:,2))
plot(t,traj_opt(:,3))
xlabel("Time")
legend("X", "Y", "\theta")
title("Ego Vehicle Pose v/s Time")

% Plotting slack variables
figure(3)

subplot(4,1,1)
plot(t(1:N), u_opt(:,3))
title('CLF Slack Variable v/s Time')

subplot(4,1,2)
plot(t(1:N), u_opt(:,4))
title('Pedestrian Slack Variable v/s Time')

subplot(4,1,3)
plot(t(1:N), u_opt(:,5))
title('Lane Keeping Slack Variable v/s Time')

subplot(4,1,4)
plot(t(1:N), u_opt(:,6))
title('Min Velocity Slack Variable v/s Time')

%% Functions

function cost = objective_fun(u)
    % Objective function for minimization
    cost_coeff_v = 0.1;
    cost_coeff_omega = 0.1;
    cost_coeff_clf_slack = 0.01;
    cost_coeff_cbf_slack_1 = 100;
    cost_coeff_cbf_slack_2 = 1.5;
    cost_coeff_cbf_slack_3 = 0.1;

    cost = cost_coeff_v*sum(u(:,1).^2) + cost_coeff_omega*sum(u(:,2).^2) + ...
           cost_coeff_clf_slack*sum(u(:,3).^2) + ...
           cost_coeff_cbf_slack_1*sum(u(:,4).^2) + ...
           cost_coeff_cbf_slack_2*sum(u(:,5).^2) + ...
           cost_coeff_cbf_slack_3*sum(u(:,6).^2);
end

function [ineq, eq] = nl_constraints(u, init_pose_e, ref_traj, ...
                        pedestrian_position,N,T,enforced_rules, dist_from_ped)
    % Non-linear constraints for optimization. Includes system dynamics,
    % CLF trajectory tracking constraints and CBF rule constraints.
    v_store = u(:,1);
    omega_store = u(:,2);
    delta_e_store = u(:,3);
    r = dist_from_ped; % radius of safe circle around pedestrian
    % Tuning parameters
    eps = 1.5; % tuning constant for CLF
    gamma1 = 0.1; % tuning constant for CBF 1
    gamma2 = 100; % tuning constant for CBF 2
    gamma3 = 1; % tuning constant for CBF 3
    
    x = zeros(N,1);
    y = zeros(N,1);
    theta = zeros(N,1);
    x(1,:) = init_pose_e(1);
    y(1,:) = init_pose_e(2);
    theta(1,:) = init_pose_e(3);
    % System dynamics using forward Euler
    for i =1:N-1
        x(i+1) = x(i) + v_store(i)*cos(theta(i))*T;
        y(i+1) = y(i) + v_store(i)*sin(theta(i))*T;
        theta(i+1) = theta(i) + omega_store(i)*T;
    end
    
    % Calculate 2-norm of error between current state trajectory and 
    % reference trajectory, which is chosen as the Lyapunov function for
    % CLF. This will be used for the CLF inequality constraint.
    error = (x - ref_traj(:,1)).^2 + (y - ref_traj(:,2)).^2 + ...
            (theta - ref_traj(:,3)).^2;
    
    % Inequality constraints for fmincon. Here are the constraints in
    % order:
    %   1. CLF constraint for trajectory tracking
    %   2. CLF slack variable > 0 for all t
    %   3. and 4. 0 m/s <= ego velocity <= 20 m/s
    %   5. and 6. -50 rad/s <= ego angular velocity <= 50 rad/s
    %   7. CBF constraint 1 - Safe distance from pedestrian
    %   8. CBF constraint 2 - Lane-keeping
    %   9. CBF constraint 3 - Min. speed limit
    %   10. Non-relaxable lane constraint
    ineq = [v_store.*(2*(x - ref_traj(:,1)).*cos(theta) + ...
            2*(y - ref_traj(:,2)).*sin(theta) + ...
            2*(theta - ref_traj(:,3)).*omega_store) + eps*(error) ...
            - delta_e_store; 
            -delta_e_store;
            -v_store;
            v_store-20;
            -omega_store-deg2rad(50);
            omega_store-deg2rad(50);
            -(x-pedestrian_position(1)).^2 - (y-pedestrian_position(2)).^2 ...
            + r^2 - ...
            gamma1*(2*(x - pedestrian_position(1)).*v_store.*cos(theta) + ...
            2*(y - pedestrian_position(2)).*v_store.*sin(theta)) + u(:,4);
            -(-v_store.*cos(theta) + gamma2*(-y + 2) - u(:,5));
            -(gamma3*(v_store - 5) - u(:,6));
            0.75-y];
    
    % Equality constraints for fmincon. Includes terminal constraint for
    % final pose of eg vehicle and enforces the condition that ego must
    % come to a stop at the end of the trajectory.
    eq = [x(N) - ref_traj(N,1);
          y(N) - ref_traj(N,2);
          theta(N) - ref_traj(N,3);
          v_store(N);
          omega_store(N)];
    
    % If rule no. 1 is enforced, enforce CBF slack variable 1 to be zero,
    % and so on for the other rules.
    if (enforced_rules(1) == 1)
        eq(end+1:end+N,1) = u(:,4);
    end
    if (enforced_rules(2) == 1)
        eq(end+1:end+N,1) = u(:,5);
    end
    if (enforced_rules(3) == 1)
        eq(end+1:end+N,1) = u(:,6);
    end
end