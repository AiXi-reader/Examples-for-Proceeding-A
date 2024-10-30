% Parameter settings
T = 1; % Time range
N = 10000; % Number of discrete time steps
dt = T/N; % Time step size
t = linspace(0, T, N); % Time vector

% Initial conditions
u0 = [2; 2]; % Initial state u
v0 = [2; 2]; % Initial condition v
z0 = [2; 2]; % Initial state z

% Continuous kappa parameters
kappa_values = linspace(0, 5000, 100); % Generate 100 kappa values between 0 and 5000
num_kappa = length(kappa_values);

% Preallocate mean square error matrix
E_squared_diff_u_14 = zeros(num_kappa, 1); % Mean square error for u and z at t=1/4 
E_squared_diff_mid_u = zeros(num_kappa, 1); %  Mean square error for u and z at t=1/2 
E_squared_diff_end_u = zeros(num_kappa, 1); %  Mean square error for u and z at t=1

% Preallocate mean square error matrix
E_squared_diff_v_14 = zeros(num_kappa, 1); % Mean square error for v and z at t=1/4 
E_squared_diff_mid_v = zeros(num_kappa, 1); % Mean square error for v and z at t=1/2
E_squared_diff_end_v = zeros(num_kappa, 1); % Mean square error for v and z at t=1u1 

% Determine indices for 1/4, 1/2, and 1
quarter_idx = round(N / 4); % 1/4 
mid_idx = round(N / 2);     % 1/2
end_idx = N;                % 1

% Record states of u, v and z
U_all = zeros(2, N); % Store the state of system u
V_all = zeros(2, N); % Store the state of system v
Z_all = zeros(2, N); % Store the state of system z

for j = 1:num_kappa
    kappa = kappa_values(j);
    
    % Preallocate solution matrices
    U = zeros(2, N); % Solution matrix for system u
    V = zeros(2, N); % Solution matrix for system v
    Z = zeros(2, N); % Solution matrix for system z
    
    % Set initial conditions
    U(:,1) = u0; 
    V(:,1) = v0;
    Z(:,1) = z0; 
    
    % Euler–Maruyama method to solve differential equations
    for i = 1:N-1
        % Noise terms
        dW1 = sqrt(dt) * randn; % Real noise dW1
        dW2 = sqrt(dt) * randn; % Real noise dW2
        
        sigma_u1 = sqrt(U(1,i)^2 + 3); % Sigma for system u1
        sigma_u2 = sqrt(U(2,i)^2 + 3); % Sigma for system u2

        sigma_v1 = sqrt(V(1,i)^2 + 4); % Sigma for system v1
        sigma_v2 = sqrt(V(2,i)^2 + 4); % Sigma for system v2
        
        % Right-hand side function for system u
        f_u1 = -U(1,i) - U(2,i) + sqrt(U(1,i)^2 + 1) - 10*U(1,i) - kappa*(U(1,i) - V(1,i));
        f_u2 = U(1,i) - U(2,i) + sqrt(U(2,i)^2 + 1) - 10*U(2,i) - kappa*(U(2,i) - V(2,i));

        % Right-hand side function for system v
        f_v1 = -V(1,i) - V(2,i) + sqrt(V(1,i)^2 + 2) - 10*V(1,i) - kappa*(V(1,i) - U(1,i));
        f_v2 = V(1,i) - V(2,i) + sqrt(V(2,i)^2 + 2) - 10*V(2,i) - kappa*(V(2,i) - U(2,i));
        
        % Right-hand side function for system z
        f_z1 = -Z(1,i) - Z(2,i) + 0.5*sqrt(Z(1,i)^2 + 1) + 0.5*sqrt(Z(1,i)^2 + 2) - 10*Z(1,i);
        f_z2 = Z(1,i) - Z(2,i) + 0.5*sqrt(Z(2,i)^2 + 1) + 0.5*sqrt(Z(2,i)^2 + 2) - 10*Z(2,i);
        
        % Calculate the next state for system u, u_1 and u_2 are perturbed by the same noise dW1
        U(:,i+1) = U(:,i) + dt * [f_u1; f_u2] + [sigma_u1; sigma_u2] * dW1;

        % Calculate the next state for system v, v_1 and v_2 are perturbed by the same noise dW2
        V(:,i+1) = V(:,i) + dt * [f_v1; f_v2] + [sigma_v1; sigma_v2] * dW2;

        % Calculate the next state for system z, z_1 and z_2 are perturbed by dW1 and dW2
        Z(:,i+1) = Z(:,i) + dt * [f_z1; f_z2] + 0.5 * [sqrt(Z(1,i)^2 + 3); sqrt(Z(2,i)^2 + 3)] * dW1 + 0.5 * [sqrt(Z(2,i)^2 + 4); sqrt(Z(2,i)^2 + 4)] * dW2;
    end

    % Store the system states for the current kappa
    U_all = U; % Store the state of system u
    V_all = V; % Store the state of system v
    Z_all = Z; % Store the state of system z

    % Calculate the Euclidean mean square error:

    % t=1/4  for u and z
    squared_diff_u_14 = sum((U(:, quarter_idx) - Z(:, quarter_idx)).^2);
    E_squared_diff_u_14(j) = squared_diff_u_14; % Save 1/4 point error

    % t=1/4 for v and z
    squared_diff_v_14 = sum((V(:, quarter_idx) - Z(:, quarter_idx)).^2);
    E_squared_diff_v_14(j) = squared_diff_v_14; % Save 1/4 point error
    
    % t=1/2  for u and z
    squared_diff_mid_u = sum((U(:, mid_idx) - Z(:, mid_idx)).^2);
    E_squared_diff_mid_u(j) = squared_diff_mid_u; % Save midpoint error

    % t=1/2  for v and z
    squared_diff_mid_v = sum((V(:, mid_idx) - Z(:, mid_idx)).^2);
    E_squared_diff_mid_v(j) = squared_diff_mid_v; % Save midpoint error
    
    % t=1  for u and z
    squared_diff_end_u = sum((U(:, end_idx) - Z(:, end_idx)).^2);
    E_squared_diff_end_u(j) = squared_diff_end_u; % Save endpoint error

    % t=1  for v and z
    squared_diff_end_v = sum((V(:, end_idx) - Z(:, end_idx)).^2);
    E_squared_diff_end_v(j) = squared_diff_end_v; % Save endpoint error
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

% Mean square error at t=1/4  for v and z
subplot(6,1,4);
plot(kappa_values, E_squared_diff_v_14, 'g-', 'LineWidth', 1.5);
xlabel('\kappa');
ylabel('Mean Squared Error');
title('Mean Squared Error at t=1/4 Point Between Subsystem v^\kappa and System z');
grid on;

% Mean square error at t=1/2  for v and z
subplot(6,1,5);
plot(kappa_values, E_squared_diff_mid_v, 'r-', 'LineWidth', 1.5);
xlabel('\kappa');
ylabel('Mean Squared Error');
title('Mean Squared Error at t=1/2 Between Subsystem v^\kappa and System z');
grid on;

% Mean square error at t=1  for v and z
subplot(6,1,6);
plot(kappa_values, E_squared_diff_end_v, 'b-', 'LineWidth', 1.5);
xlabel('\kappa');
ylabel('Mean Squared Error');
title('Mean Squared Error at t=1 Between Subsystem v^\kappa and System z');
grid on;


