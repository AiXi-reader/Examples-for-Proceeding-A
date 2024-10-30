% Paramter settings 
T = 1; % Time range
N = 10000; % Number of discrete time steps
dt = T/N; % Time step size
t = linspace(0, T, N); % Time vector

% Coefficients of system 
gamma = 3; 
a = 1; 
b = 0.5; 
g = 0.1; 
kappa_values = [0, 50, 100, 500]; % Coupling strength kappa

% Initial data
u0 = 2; % initial data for u
v0 = 2; % initial data  for v
z0 = 2; % initial data  for z

figure;
for k = 1:length(kappa_values)
    kappa = kappa_values(k); % Select the current kappa value
    
    % Preallocate solution vectors
    u = zeros(1, N); % The solution vector of system u
    v = zeros(1, N); % The solution vector of system v
    z = zeros(1, N); % The solution vector of system z

    % Set initial conditions
    u(1) = u0; 
    v(1) = v0; 
    z(1) = z0; 

    % Euler–Maruyama method to solve differential equations
    for i = 1:N-1
        % Noise terms (BM)
        dW1 = sqrt(dt) * randn; % real-valued noise dW1
        dW2 = sqrt(dt) * randn; % real-valued noise dW2
        
        % Right-hand side functions for system u and v
        du = (-gamma * u(i) - kappa * (u(i) - v(i)) - a * sin(u(i)) + g) * dt + b * sin(u(i)) * dW1;
        dv = (-gamma * v(i) - kappa * (v(i) - u(i)) - a * sin(v(i)) + g) * dt + b * sin(v(i)) * dW2;

        % Calculate the next state for system u and system v
        u(i+1) = u(i) + du;
        v(i+1) = v(i) + dv;
        
        % Right-hand side functions for system z
        dz = (-gamma * z(i) - a * sin(z(i)) + g) * dt + (b/2) * (sin(z(i)) * dW1 + sin(z(i)) * dW2);
        
        % Calculate the next state for system z
        z(i+1) = z(i) + dz;
    end

    % Plot changes of u, v and z for the current kappa
    subplot(2, 2, k); % Create 2*2 subplot
    plot(t, u, 'r', 'DisplayName', 'u^\kappa(t)');
    hold on;
    plot(t, v, 'b', 'DisplayName', 'v^\kappa(t)');
    plot(t, z, 'g', 'DisplayName', 'z(t)');
    xlabel('Time (t)');
    ylabel('Value');
    title(['Trajectories u^\kappa(t),v^\kappa(t), z(t)  for \kappa = ', num2str(kappa)]);
    legend('show');
    grid on;
end
hold off;
