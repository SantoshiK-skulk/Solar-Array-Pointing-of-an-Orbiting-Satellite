clc;
clear;
close all;

load('angle_desired')

%motor parameters
Vdc = 120;              %DC voltage 120V
a_max = 1;              %max current 1A
theta_dot_max = 0.5;    %max angular speed 0.5 deg/s
time_max = 30;          %time to reach max angular speed
t_sample = 0.5;

%derived parameters
tau = 0.6321*time_max;  %time constant: reach 63% of max speed

motor = tf([theta_dot_max],[tau 1 0]);
% figure;step(motor)
%plant estimation - this is the best I can do without knowing 
%things like inductance, resistance, etc
%http://ctms.engin.umich.edu/CTMS/index.php?aux=Activities_DCmotorB
%Plant = theta_dot_max/(tau*s + 1) * 1/s <= 1/s for theta (integration of theta_dot)
motor = c2d(motor,t_sample,'zoh');
% figure;step(motor)

num = motor.NUM{1};
den = motor.DEN{1};
% In observable form
A = [...
    -den(2) 1;...
    -den(3) 0];

B = [...
    num(2);...
    num(3)];

C = [1 0]; 

D = 0;
[a,b,c,d] = tf2ss(num,den)



plant = ss(A,B,C,D,t_sample); % in case we want a ss system
%In controllable canonical form
[an,bn,cn,dn] = tf2ss([theta_dot_max],[tau 1]);
as = [0 1;...
    0 an]; 

bs = [0;...
    bn]; cs = [cn 0]; ds=0;
plants = ss(as,bs,cs,ds);
plantsd = c2d(plants,0.5,'zoh');
[A,B,C,D] = ssdata(plantsd)
%% DESIGN PARAMETERS: Items we can change
alpha = alpha_des(:,2)'; alpha_vel = alpha_vel_des';
beta = beta_des(:,2)'; beta_vel = beta_vel_des';

% Simulation Parameters 
time = 0:t_sample:(length(alpha)*(t_sample) - t_sample); % time vector for simulation
% r = alpha;
r = beta;
% time = 0:t_sample:300; % simulate 5 minutes
% r = 10*ones(1,length(time)); % step reference signal

x0 = [(r(1)-20)/cn;0]; % True initial condition
x_hat0 = x0 + 0.05*randn(2,1); % Initial state estimate

noise_var = (0.03)^2; % noise variance

% PLACE CONTROLLER POLES
poles_desired = [0.88 0.87];
K = place(A,B,poles_desired);

G = 1/(C*inv(eye(2)-A+B*K)*B); % a.k.a. K0, for steady state gain

% LQR Control Parameters
Q = diag([1,100]); % state penalties
R = diag(1e5);     % input penalties
tspan = 0:t_sample:600;    % time span for LQR optimization
tspan = time;
%OBSERVER MODELING 
%poles of plant system ~5x faster than controller
poles = [0.1 -0.1];

%Calculate gains for observer
L = (place(A',C',poles))';

% l = place(asd',csd',poles)';
%% Simulation / Closed Loop for pole placement 

Kfun.time = 0;
Kfun.signals.values = K;
Kfun.signals.dimensions=size(K);

x = zeros(2,length(time));
x(:,1) = x0;
xhat = x;
xhat(:,1) = x_hat0;
u = zeros(1,length(time));

for idx = 1:length(time)-1
    
    u(idx) = G*r(:,idx) - K*xhat(:,idx);

    if abs(u(idx)) > 1
        u(idx) = 1*sign(u(idx));
    end
    
    n = noise_var.*randn(1,1);
    y(idx) = C*x(:,idx) + n;
    yhat(idx) = C*xhat(:,idx);
    
    xhat(:,idx+1) = A*xhat(:,idx) + B*u(idx) + L*(y(idx)-yhat(idx));
    x(:,idx+1) = A*x(:,idx) + B*u(idx);
    
end
upole = u;
ypole = [y C*x(:,idx+1)+n];
yhatpole = [yhat C*x(:,idx+1)];
figure(1); clf;

stairs(time,ypole,'LineWidth',2)
hold on
stairs(time,yhatpole,'--','LineWidth',2)
plot(time,r,'-.','LineWidth',1.5)
% axis square
title("Pole Placement Response",'FontSize',20,'Interpreter','latex')
legend('$\theta_m$','$\hat{\theta}$','$\theta_{ref}$','Interpreter','latex',...
    'Location','Southeast','FontSize',20)
xlabel("Time (s)",'FontSize',15,'Interpreter','latex')
ylabel("Angle (deg.)",'FontSize',15,'Interpreter','latex')

upole = u;
%% LQR
Afun = @(i) A;
Bfun = @(i) B;

[Kl,~] = lqr_LTV(Afun,Bfun,Q,R,tspan);

%% LQR Sim

xlqr = zeros(2,length(time));
xlqr(:,1) = x0;
xhatlqr = xlqr;
xhatlqr(:,1) = x_hat0;
u = zeros(1,length(time));

for idx = 1:length(time)-1
    idx2 = mod(idx,length(Kl)-1)+1;
    Ki = Kl{idx2};
    Gi = 1/(C*pinv(eye(2)-A+B*Ki)*B); % a.k.a. K0
    
    u(idx) = Gi*r(idx) - Ki*xhatlqr(:,idx);
    
    if abs(u(idx)) > 1
        u(idx) = 1*sign(u(idx));
    end
    
    n = noise_var.*randn(1,1);
    y(idx) = C*xlqr(:,idx) + n;
    yhat(idx) = C*xhatlqr(:,idx);

    xlqr(:,idx+1) = A*xlqr(:,idx) + B*u(idx);
    xhatlqr(:,idx+1) = A*xhatlqr(:,idx) + B*u(idx) + L*(y(idx)-yhat(idx));

end
ylqr = [y C*x(:,idx+1)+n];
yhatlqr = [yhat C*x(:,idx+1)];
ulqr = u;
figure(2); clf;
stairs(time,ylqr,'LineWidth',2)
hold on
stairs(time,yhatlqr,'--','LineWidth',2)
plot(time,r,'-.','LineWidth',1.5)
legend('$\theta_m$','$\hat{\theta}$','$\theta_{ref}$','Interpreter','latex',...
    'Location','Southeast','FontSize',20)
% axis square
title("LQR Response",'FontSize',20,'Interpreter','latex')
xlabel("Time (s)",'FontSize',15,'Interpreter','latex')
ylabel("Angle (deg.)",'FontSize',15,'Interpreter','latex')
%%
figure(3); clf;
% yyaxis right
stairs(time,upole,'LineWidth',2) 
hold on;
stairs(time,ulqr,'LineWidth',2)
% title("Pole Placement Control Input",'FontSize',20,'Interpreter','latex')
title("Control Input",'FontSize',20,'Interpreter','latex')
xlabel("Time (s)",'FontSize',15,'Interpreter','latex')
ylabel("Input",'FontSize',15,'Interpreter','latex')
legend('Pole Placement','LQR','Interpreter','latex',...
    'Location','Northeast','FontSize',20)

% figure(4); clf;
% % yyaxis right
% stairs(time,ulqr,'LineWidth',2) 
% title(" LQR Control Input",'FontSize',20,'Interpreter','latex')
% xlabel("Time (s)",'FontSize',15,'Interpreter','latex')
% ylabel("Input",'FontSize',15,'Interpreter','latex')
% % axis square

%% Compare State Errors
error_lqr = abs(ylqr - r);
error_pole_placement = abs(ypole - r);
figure(5); clf;
stairs(time,error_lqr,'LineWidth',2)
hold on;
stairs(time,error_pole_placement,'LineWidth',2)
title("Control Error",'FontSize',20,'Interpreter','latex')
xlabel("Time (s)",'FontSize',15,'Interpreter','latex')
ylabel("Error (deg.)",'FontSize',15,'Interpreter','latex')
legend('LQR Control','Pole Placement','Interpreter','latex',...
    'Location','Northeast','FontSize',20)
%% Observer Plot
time = 0:t_sample:60;
xtrue = zeros(2,length(time));
xtrue(:,1) = x0;
xobs = xtrue;
xobs(:,1) = x0+5*randn(2,1);

clear y yhat

for idx = 1:length(time)-1
    n = noise_var.*randn(1,1);
    y(idx) = C*xtrue(:,idx) + n;
    yhat(idx) = C*xobs(:,idx);

    xtrue(:,idx+1) = A*xtrue(:,idx);
    xobs(:,idx+1) = A*xobs(:,idx) + L*(y(idx)-yhat(idx));

end
y(idx+1) = C*xtrue(:,idx+1) + n;
yhat(idx+1) = C*xobs(:,idx+1);

n = noise_var.*randn(1,1);
y(idx+1) = C*xtrue(:,idx+1) + n;
figure(6); clf;
plot(time,y,'--','LineWidth',1.5);
hold on;
plot(time,yhat,'LineWidth',2);
title("Observer Verification",'FontSize',20,'Interpreter','latex')
xlabel("Time (s)",'FontSize',15,'Interpreter','latex')
ylabel("Value",'FontSize',15,'Interpreter','latex')
legend('Measured','Estimated','Interpreter','latex',...
    'Location','Northeast','FontSize',20)
%% Helper Functions
function [Acl,Bcl,Ccl,Dcl] = getClosedLoop(A,B,C,L,K,G)
Acl = [...
    A      -B*K;...
    L*C    A-B*K-L*C];

Bcl = [...
    B*G;...
    B*G];

Ccl = [C zeros(size(C))];
Dcl = zeros(size(Ccl,1),size(Bcl,2));

% Ccl_input = [-K zeros(size(C))];
% Dcl_input = G;
end