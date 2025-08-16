%% Aylin Öztürk, 20016928, Department of Control & Automation Engineering
%% PART A
clear all % Clear all variables from the workspace
clc       %Clear the command window
close all
syms  a1 a2 a3 a4 s e  % representing alphas as a symbolic value

m=8; %last digit of my student number
I=eye(4); % determine I matrix - identity matrix 4x4
%state space representation of the system: (dx/dt=A*x+B*u)
A=[0 1 0 0; 0 0 1 0; 0 0 0 1; -a4 -a3 -a2 -a1]; % matrix of system dynamics-A matrix
B=[0;0;0;1]; %input matrix-B matrix

Cont=[B A*B (A^2)*B (A^3)*B]; %our controllbility matrix
if rank(Cont) ~= 4 % checking whether system is controllable or not
    print("system is not controllable"); % printing
    error('error-not controllable'); % if it is not , throw an error
end

char_equat = (s+m).^4; % desired poles at -m, our characteristic eq.
desired_coeffs_K = coeffs(char_equat, s, 'All') %getting coefficients of char_equat
%constructing the characteristic equation matrix (phi_A) and implement A
%matrix into our phi_A
phi_A =desired_coeffs_K(1)*A^4+ desired_coeffs_K(2)*(A^3) + desired_coeffs_K(3)*(A^2) + desired_coeffs_K(4)*A + desired_coeffs_K(5)*I; 
K = [0 0 0 1] * inv(Cont) * phi_A; % calculation of K matrix- acakermann formulation for K matrix

disp('state feedback gain matrix K is:');disp(K); % result is displayed

%% PART B
A_Closed_loop= A-(B*K); %closed-loop sys. matrix- new A matrix of closed loop system
x0_states=[m/2;m/2;m;m]; %initial conditions of states are defined
simulation_time=[0:1:10]; %10saniye
x_t = zeros(4, length(simulation_time)); %  x_t matrix to store state trajectories (4 states x time points)

for i = 1:length(simulation_time) % compute state trajectories over time-10sec.
    t = simulation_time(i);%current time
    x_t(:, i) = (expm(A_Closed_loop * t)) * x0_states; %computing states using e^(A_Closed_loop * t) * x0- eA*t
end

% Plot state trajectories
figure;%create a new figure
hold on;%  hold the plot to add multiple curves
plot(simulation_time, x_t(1, :), 'r','LineWidth', 2); %x1(t)-state 1 of the closed loop system -(A-BK)
plot(simulation_time, x_t(2, :), 'g','LineWidth', 2); % x2(t)-state 2 of the closed loop system -(A-BK)
plot(simulation_time, x_t(3, :), 'b','LineWidth', 2); % x3(t)-state3 of the closed loop system -(A-BK)
plot(simulation_time, x_t(4, :), 'k','LineWidth', 2); %x4(t) -state4 of the closed loop system -(A-BK)

xlabel('time in sn');ylabel('state values'); % axis labeling
legend('x1(t)', 'x2(t)', 'x3(t)', 'x4(t)'); % adding legend to plots
title('State Trajectories of the Closed-Loop System (A-BK)'); % title of plot
grid on; %for better visualization
hold off;% elease the hold on the plot

%% PART C
C= [1 0 0 0]; %output matrix - y=C*x+D*u
Observ=[C;C*(A);C*(A^2);C*(A^3)]; % our observability matrix

if rank(Observ) ~= 4 % checking whether system is observable or not
   print("system is not observable"); % printing error
   error('error-not observable'); % if it is not , system will be throwing an error
end

p_observer = (s + 3*m).^4; % Desired characteristic polynomial for observer
desired_coeffs_observer = coeffs(p_observer, s, 'All'); %getting coefficients of p_observer
%constructing the characteristic equation matrix (phi_A) and implement A
phi_A_obs=desired_coeffs_observer(1)*A^4 +desired_coeffs_observer(2)* A^3 +desired_coeffs_observer(3)*A^2 +desired_coeffs_observer(4)*A +desired_coeffs_observer(5)*I;
L=phi_A_obs*inv(Observ)*[0;0;0;1]; % ackermann formulation for L observer gain matrix

disp('Observer gain matrix L is:');disp(L); % Display observer gain matrix L
A_observer= A-(L*C);
disp(A_observer)
%% D ŞIKKI
% implementing a1=1, a2=1, a3=1, a4=1 and creating new numeric matrises
A_new=[0 1 0 0; 0 0 1 0; 0 0 0 1; -1 -1 -1 -1]; % new numeric A system matrix
phi_A_new =A_new^4+ 32*(A_new^3) + 384*(A_new^2) + 2048*A_new+ 4096*I; %our new A matrix into our characteristic equation
Cont_new=[B A_new*B (A_new^2)*B (A_new^3)*B]; % new numeric controllability matrix
K_numeric = [0 0 0 1] * inv(Cont_new) * phi_A_new; % calculation of new K matrix- acakermann formulation for K matrix
Observ_new=[C;C*(A_new);C*(A_new^2);C*(A_new^3)]; % our observability matrix
phi_A_obs_new= desired_coeffs_observer(1)*A_new^4 + desired_coeffs_observer(2)* A_new^3  ...
+ desired_coeffs_observer(3)*A_new^2 + desired_coeffs_observer(4)*A_new +desired_coeffs_observer(5)*I;
B=[0;0;0;1];

L_new = phi_A_obs_new*inv(Observ_new)*[0;0;0;1]; % ackermann formulation for L observer gain matrix
%combining plant and observer dynamics is A_fuullmatrix:
A_fullmatrix = double([A_new-(B*K_numeric), (B*K_numeric); zeros(size(A_new)), A_new-(L_new*C)]);

x0 = [m/2; m/2; m; m]; % initial values of actual states
x0_estimated = [1; 1; 1; 1]; %initial values of estimated states
x0_full = [x0; x0_estimated]; % containing actual states and estimated states- x0_full

t = 0:0.01:10;%time vector 

x_t_full = zeros(length(x0_full), length(t)); %initialize matrix to store state evolution over time
% i am computing all state trajectories over time -10sec.
for i = 1:length(t)
    x_t_full(:, i) = expm(A_fullmatrix * t(i)) * x0_full;% Solve for states at each time 
end

x_t = x_t_full(1:4, :); %extracting actual s from full state trajectory
x_t_estimated = x_t_full(5:8, :); % extracting estimated states from full state trajectory

% x1 vs x1_estimated
figure; % create new figure for plotting
hold on; %hold the plot to add multiple curves 
plot(t, x_t(1, :), 'b','LineWidth', 3); % plot actual x1 trajectory
plot(t, x_t_estimated(1, :), '--r','LineWidth', 2);% plot estimated x1 trajectory
xlabel('time in sec.');ylabel('x1 & x1_estimated of overall system');%axis labeling
legend('x1(t)', 'x1(t)estimated'); %adding legend
title('actual and estimated state trajectories of the overall closedloop system - x1(t)'); %title of graphic
grid on;
hold off;

% x2 vs x2_estimated
figure;
hold on;
plot(t, x_t(2, :), 'b','LineWidth', 3); % plot actual x2 trajectory
plot(t, x_t_estimated(2, :), '--r','LineWidth', 2); %plot actual x2 trajectory
xlabel('time in sec.');ylabel('x2 & x2_estimated');%axis labeling
legend('x2(t)', 'x2(t)estimated of overall system');%adding legend
title('actual and estimated state trajectories of the overall closed loop system -x2(t) ');%title of graphic
grid on;hold off;

% x3 vs x3_estimated
figure;
hold on;
plot(t, x_t(3, :), 'b','LineWidth', 3);% plot actual x3 trajectory
plot(t, x_t_estimated(3, :), '--r','LineWidth', 2); %plot actual x3 trajectory
xlabel('time in sec.');ylabel('x3 & x3_estimated');%axis labeling
legend('x3(t)', 'x3(t)estimated of overall system');%adding legend
title('actual and estimated state trajectories of the overall closed loop system -x3(t)');%title of graphic
grid on;hold off;

% x4 vs x4_estimated
figure;
hold on;
plot(t, x_t(4, :), 'b','LineWidth', 3); % plot actual x4 trajectory
plot(t, x_t_estimated(4, :), '--r','LineWidth', 2); % plot actual x4 trajectory
xlabel('time in sec.');ylabel('x4 & x4_estimated of overall system'); %axis labeling
legend('x4(t)', 'x4(t)estimated');%adding legend
title('actual and estimated state trajectories of the overall closedloop system -x4(t)');%title of graphic
grid on;hold off;




%% E) Observer Design under Noise vs Kalman Filter
A_obs = A_new; 
B_obs = B; 
C_obs = C;

% Noise covariance
Q = 0.01*eye(4);  % process noise
R = 0.1;          % measurement noise

dt = 0.01;
Tend = 10;
t = 0:dt:Tend;
n = length(t);

% Initial conditions
x_true    = zeros(4,n);  x_true(:,1) = [m/2; m/2; m; m];
xhat_L    = zeros(4,n);  xhat_L(:,1) = [1;1;1;1];      % Luenberger
xhat_KF   = zeros(4,n);  xhat_KF(:,1) = [1;1;1;1];      % Kalman
P         = eye(4);
u         = zeros(1,n);


p_observer = (s + 3*m)^4;
desired_coeffs_observer = coeffs(p_observer, s, 'All');
phi_A_obs=desired_coeffs_observer(1)*A_obs^4 +desired_coeffs_observer(2)* A_obs^3 ...
    +desired_coeffs_observer(3)*A_obs^2 +desired_coeffs_observer(4)*A_obs ...
    +desired_coeffs_observer(5)*eye(4);
Observ = [C_obs; C_obs*A_obs; C_obs*A_obs^2; C_obs*A_obs^3];

rank_O = rank(Observ);

if rank_O == rank(A)
    disp('System is observable');
else
    disp(['System is NOT fully observable. Rank = ', num2str(rank_O)]);
end
L_new = phi_A_obs*inv(Observ)*[0;0;0;1];  % Luenberger gain

for k = 1:n-1
    % Real system
    dx = A_obs*x_true(:,k) + B_obs*u(k);
    x_true(:,k+1) = x_true(:,k) + dx*dt;

    % Measurement with noise
    y_meas = C_obs*x_true(:,k+1) + sqrt(R)*randn;

    % Luenberger observer update
    dxhat_L = A_obs*xhat_L(:,k) + B_obs*u(k) + L_new*(y_meas - C_obs*xhat_L(:,k));
    xhat_L(:,k+1) = xhat_L(:,k) + dxhat_L*dt;

    % Kalman filter update
    Pdot = A_obs*P + P*A_obs' - P*C_obs'*(R\C_obs)*P + Q;
    P = P + Pdot*dt;
    K_kf = P*C_obs'/R;
    dxhat_KF = A_obs*xhat_KF(:,k) + B_obs*u(k) + K_kf*(y_meas - C_obs*xhat_KF(:,k));
    xhat_KF(:,k+1) = xhat_KF(:,k) + dxhat_KF*dt;
end

% Plot comparison
state_labels = {'x1','x2','x3','x4'};
for i = 1:4
    figure; hold on; grid on;
   
    plot(t, xhat_L(i,:), '-g', 'LineWidth', 7);     % Luenberger with noise
    plot(t, x_true(i,:), 'k-', 'LineWidth', 3);       % Actual
    plot(t, xhat_KF(i,:), '--r', 'LineWidth', 2);   % Kalman
     
    xlabel('Time [s]'); ylabel(state_labels{i});
    legend('Luenberger','Actual','Kalman','Location','best','AutoUpdate','off');
    title(['Observer Performance under Noise - ', state_labels{i}]);
    hold off;
end
