%  Parameter settings
T = 1; % Time range
N = 10000; % Number of discrete time steps
dt = T/N; % Time step size
t = linspace(0, T, N); % Time vector
% Coefficients of system 
gamma = 3; 
a = 1; 
b = 0.5; 
g = 0.1; 
% Initial conditions
u0 = 2; % Initial state u
v0 = 2; % Initial state v
z0 = 2; % Initial state z

% Continuous kappa parameters
kappa_values = linspace(0, 5000, 100); % Generate 100 kappa values between 0 and 5000
num_kappa = length(kappa_values);

% Preallocate mean square error matrix
E_squared_diff_u_14 = zeros(num_kappa, 1); % Mean square error between  u and z at time t=1/4 
E_squared_diff_mid_u = zeros(num_kappa, 1);% Mean square error between  u and z at time t=1/2
E_squared_diff_end_u = zeros(num_kappa, 1);% Mean square error between  u and z at time t=1

E_squared_diff_v_14 = zeros(num_kappa, 1); % Mean square error between  v and z at time t=1/4
E_squared_diff_mid_v = zeros(num_kappa, 1);% Mean square error between  v and z at time t=1/2
E_squared_diff_end_v = zeros(num_kappa, 1);% Mean square error between  v and z at time t=1

% Determine indices for 1/4, 1/2, and 1
quarter_idx = round(N / 4); % 1/4 
mid_idx = round(N / 2);     % 1/2
end_idx = N;                % 1

for j = 1:num_kappa
    kappa = kappa_values(j);
    
    % Preallocate solution matrices
    U = zeros(1, N); % Solution matrix for system u
    V = zeros(1, N); % Solution matrix for system v
    Z = zeros(1, N); % Solution matrix for system z
    
    % Set initial condition
    U(1) = u0; 
    V(1) = v0;
    Z(1) = z0; 
    
    % Euler–Maruyama method to solve differential equations
    for i = 1:N-1
        % Noise terms
        dW1 = sqrt(dt) * randn; % Real-valued noise dW1
        dW2 = sqrt(dt) * randn; % Real-valued noise dW2
        
        % Right-hand side function for system u and system v
        du = (-gamma * U(i) - kappa * (U(i) - V(i)) - a * sin(U(i)) + g) * dt + b * sin(U(i)) * dW1;
        dv = (-gamma * V(i) - kappa * (V(i) - U(i)) - a * sin(V(i)) + g) * dt + b * sin(V(i)) * dW2;

        % Calculate the next state for system u and system v
        U(i+1) = U(i) + du;
        V(i+1) = V(i) + dv;
        
        % Right-hand side function for system z
        dz = (-gamma * Z(i) - a * sin(Z(i)) + g) * dt + (b/2) * (sin(Z(i)) * dW1 + sin(Z(i)) * dW2);
        
        % Calculate the next state for system z
        Z(i+1) = Z(i) + dz;
    end

    % Calculate the mean square error:
    E_squared_diff_u_14(j) = (U(quarter_idx) - Z(quarter_idx))^2; % t=1/4 for u and z
    E_squared_diff_mid_u(j) = (U(mid_idx) - Z(mid_idx))^2; % t=1/2 for u and z
    E_squared_diff_end_u(j) = (U(end_idx) - Z(end_idx))^2; % t=1 for u and z

    E_squared_diff_v_14(j) = (V(quarter_idx) - Z(quarter_idx))^2; % t=1/4 for v and z
    E_squared_diff_mid_v(j) = (V(mid_idx) - Z(mid_idx))^2; % t=1/2 for v and z
    E_squared_diff_end_v(j) = (V(end_idx) - Z(end_idx))^2; % t=1 for v and z
end

% Plot mean square error as a function of kappa
figure;

% Mean square error at t=1/4 for u and z
subplot(6,1,1);
plot(kappa_values, E_squared_diff_u_14, 'g-', 'LineWidth', 1.5);
xlabel('\kappa');
ylabel('Mean Squared Error');
title('Mean Squared Error at t=1/4 Point Between Subsystem u^\kappa and System z');
grid on;

% Mean square error at t=1/2 for u and z
subplot(6,1,2);
plot(kappa_values, E_squared_diff_mid_u, 'r-', 'LineWidth', 1.5);
xlabel('\kappa');
ylabel('Mean Squared Error');
title('Mean Squared Error at t=1/2 Between Subsystem u^\kappa and System z');
grid on;

% Mean square error at t=1 for u and z
subplot(6,1,3);
plot(kappa_values, E_squared_diff_end_u, 'b-', 'LineWidth', 1.5);
xlabel('\kappa');
ylabel('Mean Squared Error');
title('Mean Squared Error at t=1 Between Subsystem u^\kappa and System z');
grid on;

% Mean square error at t=1/4 for v and z
subplot(6,1,4);
plot(kappa_values, E_squared_diff_v_14, 'g-', 'LineWidth', 1.5);
xlabel('\kappa');
ylabel('Mean Squared Error');
title('Mean Squared Error at t=1/4 Point Between Subsystem v^\kappa and System z');
grid on;

% Mean square error at t=1/2 for v and z
subplot(6,1,5);
plot(kappa_values, E_squared_diff_mid_v, 'r-', 'LineWidth', 1.5);
xlabel('\kappa');
ylabel('Mean Squared Error');
title('Mean Squared Error at t=1/2 Between Subsystem v^\kappa and System z');
grid on;

% Mean square error at t=1 for v and z
subplot(6,1,6);
plot(kappa_values, E_squared_diff_end_v, 'b-', 'LineWidth', 1.5);
xlabel('\kappa');
ylabel('Mean Squared Error');
title('Mean Squared Error at t=1 Between Subsystem v^\kappa and System z');
grid on;
