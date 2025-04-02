%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%          Magnetic friction of rotor arrays                              %
%                                                                         %
%          Author:   Anton Lueders                                        %
%          Date:     03/2024                                              %
%          Modified by Hongri Gu                                          %
%                                                                         %
%          AG Nielaba                                                     %
%          Statistical and computational physics                          %
%          Univesity of Konstanz                                          %
%                                                                         %
%          AG Franosch                                                    %
%          BioNano-Physics                                                %
%          Univesity of Innsbruck                                         %
%                                                                         %
%          Last edited: 05.03.2025                                        %
%                                                                         %
%          V2.0.0                                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

set(groot,'DefaultTextInterpreter' ,'LaTeX');
set(groot,'DefaultAxesTickLabelInterpreter' ,'LaTeX');
set(groot,'DefaultAxesFontName' ,'LaTeX');
set(groot,'DefaultLegendInterpreter' ,'LaTeX');

rng(33);

folderPath=pwd;

%%=======================================================================%%
% Input parameter
%========================================================================%%


% Micromagnet specific parameters ========================================%

friction_coefficient = 2.5*10^(-6); % [kg*m^2/s]
moment_of_inertia = 7.48*10^(-10);  % [kg*m^2]
magnetic_moment = 4.57*10^(-2); % [A*m^2]

% Probe parameters =======================================================%

dim_x = 7;              % Number of spin in x-direction
dim_y = 7;              % Number of spin in y-direction
lattice_constant_x = 0.016; % [m] Distance between spins in x-direction
lattice_constant_y = 0.016; % [m] Distance between spins in y-direction
v_x = 0.008;              % [m/s] Translational velocity in x-direction
v_y = 0.000;              % [m/s] Translational velocity in y-direction

% Substrat parameters ====================================================%

substrat_dim_x = 4;     % Size in x-direction (multiples of dim_x)
substrat_dim_y = 4;     % Size in y-direction (multiples of dim_y)
angle_substrat = 0;     % Orientation of micromagnets 
magnetic_moment_sub = 5.03*10^(-2); % [Am^2] magnetic moment of substrat 
tilt_substrat = 0.0;

% Relative parameters ====================================================%

shift_x = 0;            % Shift in x_direction relative to the substrat
shift_y = 0.000;        % [m] Shift in y_direction relative to the substrat
shift_z = 0.0075;       % [m] Shift in z_direction relative to the substrat
tilt_angle = 0.0;       % Angle between substrat and probe

% Simulation parameters ==================================================%

time_step = 0.0003;     % [s] length of simulation step
max_time = 10.5;        % [s] Total simulation time
start_t_time = 0.5;     % [s] Time before the translation starts 

% Calc Force =============================================================%

calc_force = 1;            % If = 1, The force is calculated and printed
Force_calc_rate = 0.001;   % [s] Rate for printing force in files

scan_value=0.0050:0.0001:0.013; % Varied paramter of the simulations

for scan=1:length(scan_value)
    
    shift_z=scan_value(scan);
    
%%=======================================================================%%
% Init of variables
%%=======================================================================%%
    
    magnet_array_0 = zeros(dim_x, dim_y);
    magnet_array_0(1:2:dim_x,:) = 0 + ...
        pi/24 * (rand(length(1:2:dim_x),dim_y) - 0.5);
    magnet_array_0(2:2:dim_x,:) = pi +...
        pi/24 * (rand(length(2:2:dim_x),dim_y) - 0.5);
    velocity_array_0 = zeros(dim_x, dim_y);   % Initial velocity = 0
    torques_0 = zeros(dim_x, dim_y);
    
    magnet_array_1 = magnet_array_0 + velocity_array_0 * time_step ...
        + 1 / 2 * torques_0 * time_step^2;    % Needed for Verlet algorithm
    
    angular_velocity = CalcVelocity(magnet_array_0, ... 
            magnet_array_1, time_step);           
    
    substrat_array = angle_substrat * ones((substrat_dim_y + 1) * dim_x,...
        (substrat_dim_x + 1) * dim_y);
    
% Set coordinates of the spins on the probe ==============================%
    
    pos_in_grid_x = zeros(dim_x, dim_y);
    for x = 1:dim_x
        pos_in_grid_x(:,x) = (x - 1) * lattice_constant_x + shift_x;
    end
    
    pos_in_grid_y = zeros(dim_x, dim_y);
    for y = 1:dim_y
        pos_in_grid_y(y,:) = - (y - 1) * lattice_constant_y ...
            + (dim_y - 1) * lattice_constant_y + shift_y;
    end
    
    pos_in_grid_x_old = pos_in_grid_x;
    pos_in_grid_x = pos_in_grid_x * cos(tilt_angle) ...
        - pos_in_grid_y * sin(tilt_angle);         % Applies the tilt angle
    pos_in_grid_y = pos_in_grid_x_old * sin(tilt_angle) ...
        + pos_in_grid_y * cos(tilt_angle);         % Applies the tilt angle
    
% Set coordinates of the spins on the substrat ===========================%
    
    pos_substrat_x = zeros((substrat_dim_y + 1) * dim_y, ...
        (substrat_dim_x + 1) * dim_x);
    for x = 1:(substrat_dim_x + 1) * dim_x
        pos_substrat_x(:,x) = (x - 1) * lattice_constant_x - ...
            (substrat_dim_x / 2) * dim_x * lattice_constant_x;
    end
    
    pos_substrat_y = zeros((substrat_dim_y + 1) * dim_x, ...
        (substrat_dim_x + 1) * dim_y);
    for y = 1:(substrat_dim_y + 1) * dim_y
        pos_substrat_y(y,:) = - (y - 1) * lattice_constant_y + ...
            ((substrat_dim_y / 2) * dim_y + ...
            (dim_y - 1)) * lattice_constant_y;
    end
    
    pos_substrat_x_old = pos_substrat_x;
    pos_substrat_x = pos_substrat_x * cos(tilt_substrat) ...
        - pos_substrat_y * sin(tilt_substrat);  
    pos_substrat_y = pos_substrat_x_old * sin(tilt_substrat) ...
        + pos_substrat_y * cos(tilt_substrat); 
    
    
% Variables for saving temporary results =================================%
    
    magnet_array_0_temp = zeros(dim_x, dim_y);
    magnet_array_1_temp = zeros(dim_x, dim_y);
    pos_in_grid_x_temp = zeros(dim_x, dim_y);
    pos_in_grid_y_temp = zeros(dim_x, dim_y);
    
    data = {};     % Initialize a cell array to store the data
    
%%=======================================================================%%
% Main simulation loop
%%=======================================================================%%
      
    tic            % For calculation of the simulation time
    
    for i = time_step:time_step:max_time    
    
% Calculates the torques =================================================%
    
        torques = CalcTorque(magnet_array_1, angular_velocity, ...
            friction_coefficient, pos_in_grid_x, pos_in_grid_y, ...
            substrat_array, pos_substrat_x, pos_substrat_y, shift_z, ...
            tilt_angle, magnetic_moment, magnetic_moment_sub);
    
% Integrates the equations of motion =====================================% 
    
        [magnet_array_0_temp, magnet_array_1_temp] = ...
            Move(magnet_array_0,magnet_array_1, torques, time_step, ...
            moment_of_inertia);
    
        magnet_array_0 = magnet_array_0_temp;
        magnet_array_1 = magnet_array_1_temp;
    
% Approximation of the velocity ==========================================%
    
        angular_velocity = CalcVelocity(magnet_array_0, ... 
            magnet_array_1, time_step);
    
% Translation of the probe ===============================================%
    
        if i > start_t_time
    
            [pos_in_grid_x_temp, pos_in_grid_y_temp] ...
                = Translate(pos_in_grid_x, pos_in_grid_y, v_x, v_y, ...
                time_step);
            pos_in_grid_x = pos_in_grid_x_temp;
            pos_in_grid_y = pos_in_grid_y_temp;
    
        end
    
% Calculates the force and other observables =============================%
    
        if calc_force == 1        
            if mod(round((i+time_step) / time_step), ...  
                    round(Force_calc_rate / time_step)) == 0
                
                % Magnetic interaction force between substrate and slider
                [Fx, Fy, Fz] = CalcForce(magnet_array_1, ...
                    substrat_array, tilt_angle, magnetic_moment, ...
                    magnetic_moment_sub, shift_z, pos_in_grid_x, ...
                    pos_in_grid_y, pos_substrat_x, pos_substrat_y);
    
                % Torque acting on each moment
                T = CalcTorqueGrid(magnet_array_1, substrat_array, ...
                    tilt_angle, magnetic_moment, magnetic_moment_sub, ...
                    shift_z, pos_in_grid_x, pos_in_grid_y, ...
                    pos_substrat_x, pos_substrat_y);
    
                % Information to compute the hysteresis curves
                % Split between two sublattices phi and theta
                [T_phi, T_theta, Bx, By, Bz, Mx_phi, My_phi, Mz_phi, ...
                    Mx_theta, My_theta, Mz_theta] ...
                    = CalcHysteresis(magnet_array_1, substrat_array, ...
                    tilt_angle, magnetic_moment, magnetic_moment_sub, ...
                    shift_z, pos_in_grid_x, pos_in_grid_y, ...
                    pos_substrat_x, pos_substrat_y);
    
                % Magentic energy of the internal interactions in the
                % slider and the magnetic interaction energy between slider
                % and substrate
                [E_sub, E_probe] = CalcEnergy(magnet_array_1, ...
                    substrat_array, tilt_angle, magnetic_moment, ...
                    magnetic_moment_sub, shift_z, pos_in_grid_x, ...
                    pos_in_grid_y, pos_substrat_x, pos_substrat_y);
           end
        end
    
% Saves data =============================================================%

        if mod(round((i+time_step)/ time_step), ... 
                   round(Force_calc_rate / time_step)) == 0
	           data{end + 1} = struct('height', scan_value(scan),...
                   'time', i + time_step,...
                   'array', wrapTo2Pi(magnet_array_1),...
                   'force',[Fx, Fy, Fz],...
                   'torque_grid',T,...
                   'hysteresis',[T_phi, T_theta, Bx, By, Bz,... 
                       Mx_phi, My_phi, Mz_phi, Mx_theta, My_theta, ...
                       Mz_theta], ...
                   'energy',[E_probe, E_sub]);
        end
    end
    
    toc            % For calculation of simulation time

    % Writes save file:
    matFileName = fullfile(folderPath, ...
        ['height',num2str(scan_value(scan), '%.4f '), '.mat']);
    save(matFileName, 'data');

end