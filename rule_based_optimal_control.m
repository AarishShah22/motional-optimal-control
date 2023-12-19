clc
clear all
close all

%% Define Parameters

a = 1; % rear axle to COG (centre of gravity)
L = 2; % total length of vehicle
tf = 2; % simulation end time, simulation begins at t=0
T = 0.05; % time step
N = tf/T; % no. of time steps
x_e = 1;
y_e = 2.2;
theta_e = 0;
% init_pose_e = [1,2.2,0]; % ego vehicle initial pose (pose = [x; y; theta]) 
% ref_traj = [linspace(1,20,N)', 1.75*ones(N,1), zeros(N,1)];
x_r = linspace(1,20,N)';
y_r = 1.75*ones(N,1);
theta_r = zeros(N,1);
Xr = [x_r; y_r; theta_r];
Xp = [1.5,2.75];
priority_structure = [1,1,1;
                      1,1,0;
                      1,0,1;
                      1,0,0;
                      0,1,1;
                      0,1,0;
                      0,0,1;
                      0,0,0;];

for i = 1:8
    i
    rules = priority_structure(i,:);
    u0 = [10*ones(N,1) ones(N,1)*0.3 ones(N,1)*10 ones(N,1)*0 ones(N,1)*1 ones(N,1)*0];
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp','StepTolerance',1e-6, 'MaxFunctionEvaluations',1e6, 'MaxIterations',1e3);
    [u_opt, fval, exit_flag] = fmincon(@(u) objective_fun(u, [x_e,y_e,theta_e], Xr,N,T), u0, [], [], [], [], [], [], @(u) nl_constraints(u, [x_e,y_e,theta_e], [x_r,y_r,theta_r], Xp,N,T,rules), options);
    if exit_flag == 1 || exit_flag == 2
        break;
    end
end

%%
x(1,:) = x_e;
y(1,:) = y_e;
theta(1,:) = theta_e;

for i =1:N
    x(i+1) = x(i) + u_opt(i,1)*cos(theta(i))*T;
    y(i+1) = y(i) + u_opt(i,1)*sin(theta(i))*T;
    theta(i+1) = theta(i) + u_opt(i,2)*T;
end

t = linspace(0,tf,N+1);

figure()
plot(x,y)
hold on
scatter(x,y)
axis square
axis equal
xlim([0, 20]); % Set x-axis limits
ylim([0, 6]); % Set y-axis limits
plot([0,20],[2,2],'Linewidth',2)
viscircles([1.5,2.75],[1])

title('Trajectory Generated')

figure
plot(t,x)
hold on
plot(t,y)

figure
% plot(t(1:N), u_opt(:,1))
% hold on
% plot(t(1:N), u_opt(:,2))
% plot(t(1:N), u_opt(:,3))
plot(t(1:N), u_opt(:,5),'Linewidth',2)
title('Lane Keeping Slack Variable')

figure
plot(t(1:N), u_opt(:,4),'Linewidth',2)
title('Pedestrian Slack Variable')

figure
plot(t(1:N), u_opt(:,6),'Linewidth',2)
title('Min Velocity Slack Variable')

figure
plot(t(1:N), u_opt(:,3),'Linewidth',2)
title('CLF Slack Variable')

% figure
% plot(x,y)
% hold on
% scatter(x,y)
% plot([1,20],[2,2])
% viscircles([1.75,2],[1])
% 
% figure
% plot(t,x)
% hold on
% plot(t,y)
% 
% figure
% plot(t(1:N), u_opt(:,1))
% hold on
% plot(t(1:N), u_opt(:,2))
% plot(t(1:N), u_opt(:,3))
% legend("v","w","\delta_e")

% u0 = u_opt;

%% Functions

function cost = objective_fun(u,X0, Xr,N,T)
    coeffV = 0.1;
    coeffO = 0.1;
    coeffD = 0.01;
    coeffD1 = 100;
    coeffD2 = 1.5;
    coeffD3 = 0.1;
%     coeffError = 1;
%     x(1,:) = X0(1);
%     y(1,:) = X0(2);
%     theta(1,:) = X0(3);
%     v_store = u(:,1);
% %     size(v_store)
%     omega_store = u(:,2);
% 
%     for i =1:N
%         x(i+1) = x(i) + v_store(i)*cos(omega_store(i))*T;
%         y(i+1) = y(i) + v_store(i)*sin(omega_store(i))*T;
%         theta(i+1) = theta(i) + omega_store(i)*T;
%     end
% 
%     error = sqrt((x - Xr(1)).^2 + (y - Xr(2)).^2 + (theta - Xr(3)).^2);
    cost = coeffV*sum(u(:,1).^2) + coeffO*sum(u(:,2).^2) + coeffD*sum(u(:,3).^2) + coeffD1*sum(u(:,4).^2) + coeffD2*sum(u(:,5).^2) + coeffD3*sum(u(:,6).^2);

end

function [ineq, eq] = nl_constraints(u, X0, Xr, Xp,N,T, rules)
    
    v_store = u(:,1);
%     size(v_store)
    omega_store = u(:,2);
    delta_e_store = u(:,3);
    esp = 1.5;

    gamma1 = 0.1;
    gamma2 = 100;
    gamma3 = 1;
    x = zeros(N+1,1);
    y = zeros(N+1,1);
    theta = zeros(N+1,1);
    r = 1;

    x(1,:) = X0(1);
    y(1,:) = X0(2);
    theta(1,:) = X0(3);

    for i =1:N
        x(i+1) = x(i) + v_store(i)*cos(theta(i))*T;
        y(i+1) = y(i) + v_store(i)*sin(theta(i))*T;
        theta(i+1) = theta(i) + omega_store(i)*T;
    end
    
    x = x(1:N);
    y = y(1:N);
    theta = theta(1:N);
    error = (x - Xr(:,1)).^2 + (y - Xr(:,2)).^2 + (theta - Xr(:,3)).^2;

    ineq = [v_store.*(2*(x - Xr(:,1)).*cos(theta) + 2*(y - Xr(:,2)).*sin(theta) + 2*(theta - Xr(:,3)).*omega_store) + esp*(error) - delta_e_store ; 
        -delta_e_store;
        -v_store;
        v_store-20;
        -omega_store-deg2rad(50);
        omega_store-deg2rad(50);
        -(x-Xp(1)).^2 - (y-Xp(2)).^2 + r^2 - gamma1*(2*(x - Xp(1)).*v_store.*cos(theta) + 2*(y - Xp(2)).*v_store.*sin(theta)) + u(:,4);
        -(-v_store.*cos(theta) + gamma2*(-y + 2) - u(:,5));
        -(gamma3*(v_store - 5) - u(:,6));
        0.75-y;
        ];

    eq = [x(N) - Xr(N,1);
        y(N) - Xr(N,2);
        theta(N) - Xr(N,3);
%         v_store(N);
        omega_store(N);
        ];

    if (rules(1) == 1)
        eq(end+1:end+N,1) = u(:,4);
    end
    if (rules(2) == 1)
        eq(end+1:end+N,1) = u(:,5);
    end
    if (rules(3) == 1)
        eq(end+1:end+N,1) = u(:,6);
    end

%     eq = [x(N) - Xr(1);
%         y(N) - Xr(2);
%         theta(N) - Xr(3);
% %         v_store(N);
%         omega_store(N);
% %         u(:,3);
% %         u(:,4);
% %         u(:,5);
% %         u(:,6)
%         ];
%     if(idx == 1) % No relaxation
%         eq = [x(N) - Xr(N,1);
%         y(N) - Xr(N,2);
%         theta(N) - Xr(N,3);
%         v_store(N);
%         omega_store(N);
%         u(:,4);
%         u(:,5);
%         u(:,6)
%         ];
%     end
%     
%     if(idx == 2) % Relaxing 3
%         eq = [x(N) - Xr(N,1);
%         y(N) - Xr(N,2);
%         theta(N) - Xr(N,3);
%         v_store(N);
%         omega_store(N);
%         u(:,4);
%         u(:,5)
%         ];
%     end
%     
%     if(idx == 3) % Relaxing 2
%         eq = [x(N) - Xr(N,1);
%         y(N) - Xr(N,2);
%         theta(N) - Xr(N,3);
%         v_store(N);
%         omega_store(N);
%         u(:,4);
%         u(:,6)
%         ];
%     end
%     
%     if(idx == 4) % Relaxing 2 and 3
%         eq = [x(N) - Xr(N,1);
%         y(N) - Xr(N,2);
%         theta(N) - Xr(N,3);
%         v_store(N);
%         omega_store(N);
%         u(:,4)
%         ];
%     end
%     
%     if(idx == 5) % Relaxing 1
%         eq = [x(N) - Xr(N,1);
%         y(N) - Xr(N,2);
%         theta(N) - Xr(N,3);
% %         v_store(N);
%         omega_store(N);
%         u(:,5)
%         u(:,6)
%         ];
%     end
%     
%     if(idx == 6) % Relaxing 1 and 3
%         eq = [x(N) - Xr(N,1);
%         y(N) - Xr(N,2);
%         theta(N) - Xr(N,3);
% %         v_store(N);
%         omega_store(N);
%         u(:,5);
%         ];
%     end
%     
%     if(idx == 7) % Relaxing 2 and 1
%          eq = [x(N) - Xr(N,1);
%         y(N) - Xr(N,2);
%         theta(N) - Xr(N,3);
% %         v_store(N);
%         omega_store(N);
%         u(:,6);
%         ];
%     end
%     
%     if(idx == 8) % All relaxed
%         eq = [x(N) - Xr(N,1);
%         y(N) - Xr(N,2);
%         theta(N) - Xr(N,3);
% %         v_store(N);
%         omega_store(N);
%         ];
%     end
end