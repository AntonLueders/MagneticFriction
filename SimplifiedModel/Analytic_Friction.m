close all;
clear all;

set(groot,'DefaultTextInterpreter' ,'LaTeX');
set(groot,'DefaultAxesTickLabelInterpreter' ,'LaTeX');
set(groot,'DefaultAxesFontName' ,'LaTeX');
set(groot,'DefaultLegendInterpreter' ,'LaTeX');

diss_energy = [];
order = [];
E_Stokes = [];

% Different heights (in m) for which the simplified system is solved.
c_int = 0.0065:0.000025:0.0115;

for k = 1:length(c_int)
    
    % Time for which the system is solved and initial conditions
    tspan = 0:0.0001:10;
    angle0 = pi * [1/4, 3/4];

    % Magnitudes (in microjoule) of the three terms of 
    % the simplified Hamiltonian
    % Values are matched to the experimental setup using the relations
    % derived in the supporting information 
    a = 3 * 21 * 50.99;
    b = 42 * 50.99 + 18/sqrt(2) * 50.99;
    c = 0.5 * 42 * 2.2987 * 10^(-4) * 1/c_int(k)^3;
    
    % Angular velocity omega = k*v, matched to the experiments
    omega = pi;
    % Heuristic friction coefficient of the shaft friction. Value is chosen
    % for maximum stability and to see the three regimes
    gamma = 4.7;
    
    ode = @(t,angle_data) angle_diff_eq(t, angle_data, a, b, c, ...
        omega, gamma);
    [t, angle_data] = ode45(ode, tspan, angle0);
    
    % Rolling average to remove "numerical noise". Can be left out-
    angle_data = movmean(angle_data,1001);
    
    angle_data = angle_data(t>4*pi/omega & t <= 6*pi/omega,:);
    t = t(t>4*pi/omega & t <= 6*pi/omega);
    
    % Calculation of the substrate torque
    friction_phi = c * sin(angle_data(:,1) + omega * (t-2*pi/omega));
    friction_theta = c * sin(angle_data(:,2) + omega * (t-2*pi/omega));

    % Calculation of the order parameter
    order_parameter = sin(angle_data(:,1)) .* sin(angle_data(:,2)) ...
        + cos(angle_data(:,1)) .* cos(angle_data(:,2));
    avg_order_parameter = sum(order_parameter) / length(order_parameter);

    Delta_Angle = angle_data(2:end,:) - angle_data(1:end-1,:);
    % Rolling average to remove "numerical noise". Can be left out.
    Delta_Angle = movmean(Delta_Angle, 501);
 
    % Calculation of the dissipated energy
    int_phi = 0;
    int_theta = 0;
    for i = 1:(length(t)-1)
        int_phi = int_phi + friction_phi(i) * Delta_Angle(i,1);
        int_theta = int_theta + friction_theta(i) * Delta_Angle(i,2);
    end

    % To check the validity of the results, the friction is also
    % approximated via the Stokes law: Friction = - Gamma * omega (here, 
    % omega is called v); E_diss = - Gamma int omega^2 dt. 
    v = [];
    for i = 1:length(t)
        v = [v; angle_diff_eq(t(i), angle_data(i,:),...
            a, b, c, omega, gamma)'];
    end
    E_Stokes_phi = -gamma*sum(v(1:end-1,1).*Delta_Angle(:,1));
    E_Stokes_theta = -gamma*sum(v(1:end-1,2).*Delta_Angle(:,2));

    % Saving the data for the different heights
    E_Stokes = [E_Stokes, E_Stokes_phi + E_Stokes_theta];    
    diss_energy(k) = int_phi + int_theta;
    order(k) = avg_order_parameter;

end

% Plot of the magnetic friction calculated from the simplified model
% The line is the result of the definition using the substrate torque
% The points are the results using Stokes law.
figure(1)
plot(c_int * 10^3 ,-diss_energy * 10^(-6) / 0.016,...
    'Color',[77 195 255] / 255,'LineWidth', 2.5)
hold on;
plot(c_int * 10^3 ,-E_Stokes * 10^(-6) / 0.016,'o',...
    'Color',[77 195 255] / 255)
set(gcf, 'position',[0 0 600 400]);
set(gca, 'Layer', 'top');
set(gca,'FontSize', 21);
ax = gca;
ax.YRuler.TickLabelFormat = '%.1f';
ax.XAxis.LineWidth = 2;
ax.XAxis.Color = 'k';
ax.YAxis.Color = 'k';
ax.YAxis.LineWidth = 2;
xlabel('$ h $ [mm]','interpreter','latex', 'FontSize', 21);
ylabel('$ -<F_x^{\textrm{mag}}> $ [N]',...
    'interpreter','latex', 'FontSize', 21);
print -dpng SimplifiedModelResult;

% Saves the dissipated energy
writematrix(num2str([c_int' * 10^3, -diss_energy' * 10^(-6) / 0.016],...
    '%.4f '),'AvgForceSimplifiedModel.dat',...
         'Delimiter', 'tab')

% Saves the average order parameter
writematrix(num2str([c_int' * 10^3, order'],...
    '%.4f '),'ycoupleSimplifiedModel.dat',...
         'Delimiter', 'tab')