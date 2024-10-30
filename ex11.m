 % Parameter settings
T = 1; % Time range
N = 10000; % Number of discrete time steps
dt = T/N; % Time step size
t = linspace(0, T, N); % Time vetor

% Initial conditions 
u0 = [2; 2]; % Initial data  u
v0 = [2; 2]; % Initial data  v
z0 = [2; 2]; % Initial data  z


% kappa values to test
kappa_values = [0, 50, 500, 1000]; % Specific kappa values
num_kappa = length(kappa_values);

% Record states of u, v and z
V_all = zeros(2, N, num_kappa); % Store states of system v for each kappa
U_all = zeros(2, N, num_kappa); % Store states of system u for each kappa
Z_all = zeros(2, N, num_kappa); % Store states of system z for each kappa
for j = 1:num_kappa
    kappa = kappa_values(j);
    
    % Preallocate solution matrices
    V = zeros(2, N);% Solution matrix for system v
    U = zeros(2, N);% Solution matrix for system u
    Z = zeros(2, N);% Solution matrix for system z
    
    % Set initial data
    V(:,1) = v0;
    U(:,1) = u0; 
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
        

         % Right-hand side functions for system u
        f_u1 = -U(1,i) - U(2,i) + sqrt(U(1,i)^2 + 1) - 10*U(1,i) - kappa*(U(1,i) - V(1,i));
        f_u2 = U(1,i) - U(2,i) + sqrt(U(2,i)^2 + 1) - 10*U(2,i) - kappa*(U(2,i) - V(2,i)); 

        % Right-hand side functions for system v
        f_v1 = -V(1,i) - V(2,i) + sqrt(V(1,i)^2 + 2) - 10*V(1,i) - kappa*(V(1,i) - U(1,i));
        f_v2 = V(1,i) - V(2,i) + sqrt(V(2,i)^2 + 2) - 10*V(2,i) - kappa*(V(2,i) - U(2,i));
        
        % Right-hand side functions for system z
        f_z1 = -Z(1,i) - Z(2,i) + 0.5*sqrt(Z(1,i)^2 + 1) + 0.5*sqrt(Z(1,i)^2 + 2) - 10*Z(1,i);
        f_z2 = Z(1,i) - Z(2,i) + 0.5*sqrt(Z(2,i)^2 + 1) + 0.5*sqrt(Z(2,i)^2 + 2) - 10*Z(2,i);
        
        % Calculate the next state of system u, u_1 and u_2 are disturbed by the same noise dW1
        U(:,i+1) = U(:,i) + dt * [f_u1; f_u2] + [sigma_u1; sigma_u2] * dW1;
        % Calculate the next state of system v, v_1 and v_2 are disturbed by
        % the same noise dW2
        V(:,i+1) = V(:,i) + dt * [f_v1; f_v2] + [sigma_v1; sigma_v2] * dW2;

        % Calculate the next state of system z, z_1 and z_2 are disturbed by dW1 and dW2 
        Z(:,i+1) = Z(:,i) + dt * [f_z1; f_z2] + 0.5 * [sqrt(Z(1,i)^2 + 3); sqrt(Z(2,i)^2 + 3)] * dW1 + 0.5 * [sqrt(Z(1,i)^2 + 4); sqrt(Z(2,i)^2 + 4)] * dW2;
    end

    %  Store the current state of the system for the current kappa
    U_all(:,:,j) = U; % Store the state of system u
    V_all(:,:,j) = V; % Store the state of system v
    Z_all(:,:,j) = Z; % Store the state of system z
end

%  Plot changes of u_1 and z_1
figure;
for j = 1:num_kappa
    subplot(2, 2, j); % Create subplot
    hold on;
    plot(t, U_all(1,:,j), 'r', 'DisplayName', ['u_1^\kappa ' ]);
    plot(t, Z_all(1,:,j), 'b', 'DisplayName', ['z_1']);
    xlabel('Time (t)');
    ylabel('Value');
    title(['u_1^\kappa(t) and z_1(t) for \kappa = ', num2str(kappa_values(j))]);
    legend('show');
    grid on;
end
hold off;

%  Plot changes of u_2 and z_2
figure;
for j = 1:num_kappa
    subplot(2, 2, j); % Create subplot
    hold on;
    plot(t, U_all(2,:,j), 'g', 'DisplayName', ['u_2^\kappa']);
    plot(t, Z_all(2,:,j), 'm', 'DisplayName', ['z_2']);
    xlabel('Time (t)');
    ylabel('Value');
    title(['u^\kappa_2(t) and z_2(t) for \kappa = ', num2str(kappa_values(j))]);
    legend('show');
    grid on;
end
% Plot changes of v_1 and z_1
figure;
for j = 1:num_kappa
    subplot(2, 2, j); % Create subplot
    hold on;
    plot(t, V_all(1,:,j), 'r', 'DisplayName', ['v_1^\kappa ' ]);
    plot(t, Z_all(1,:,j), 'b', 'DisplayName', ['z_1']);
    xlabel('Time (t)');
    ylabel('Value');
    title(['v_1^\kappa(t) and z_1(t) for \kappa = ', num2str(kappa_values(j))]);
    legend('show');
    grid on;
end
hold off;

% Plot changes of v_2 and z_2
figure;
for j = 1:num_kappa
    subplot(2, 2, j); % Create subplot
    hold on;
    plot(t, V_all(2,:,j), 'g', 'DisplayName', ['v_2^\kappa']);
    plot(t, Z_all(2,:,j), 'm', 'DisplayName', ['z_2']);
    xlabel('Time (t)');
    ylabel('Value');
    title(['v^\kappa_2(t) and z_2(t) for \kappa = ', num2str(kappa_values(j))]);
    legend('show');
    grid on;
end
hold off;
