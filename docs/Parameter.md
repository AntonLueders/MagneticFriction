As written in the [Readme](../README.md), the input parameters of the simulations are written within the header of the main file **MagneticArray.m** and listed in the following form:

```
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

```

.
.
.

```
scan_value = 0.0050:0.0001:0.013; % Varied paramter of the simulations

```
.
.
.

All parameter presented in the example above are the ones utilized to perform the simulations discussed in [1]. Please check the corresponding main text and the supporting information of [1] for further information. 

The particular parameters that must be set to succesfully perform a simulation with **MagneticFriction** have the following meaning:

- **friction_coefficient:** This is the friction coefficient of the Stokes-like shaft friction. It is a free parameter of the simulations and should be chosen to maintain stability, match the experimental friction, and sustain overdamped dynamics.
- **moment_of_inertia:** This parameter is the moment of inertia of the magnetic moments of the slider. The value used in [1] is estimated by assuming a homogeneous mass and utilizing the ring magnet geometry of the experimental setup.
- **magnetic_moment:** Magnitude of the magnetic moments on the slider. For [1], this value is approximated by using the magnetization of the material of the magnets used in the experiments and the volume of the corresponding ring geometry.
- **dim_x:** Number of magnetic moments on the slider in the x-direction. For all simulations in [1], *dim_x* and *dim_y* are identical, and the source code is only tested for this setup.
- **dim_y:** Number of magnetic moments on the slider in the y-direction. For all simulations in [1], *dim_x* and *dim_y* are identical, and the source code is only tested for this setup.
- **lattice_constant_x:** Distance between magnetic moments on the slider in the x-direction.
- **lattice_constant_y:** Distance between magnetic moments on the slider in the y-direction. The values are matched to the experiments. The code is only tested for square lattices.
- **v_x:** Translational velocity in x-direction. This value is matched to the sliding speed of the experiments in [1].
- **v_y:** Translational velocity in y-direction. Should be set to zero to reproduce the results of [1]
- **substrat_dim_x:** Number of magnetic moments on the substrate in multiples of *dim_x*. Always use an even value (the code is not tested for odd values). The exact number of substrate moments in x-direction can be calculated by (*substrate_dim_x* + 1) x *dim_x*. Note that the source code is only tested for systems where *substrat_dim_x* is equal to *substrat_dim_y*. The number of magnets on the substrate is not matched to the experiments but freely chosen.
- **substrat_dim_y:** Number of magnetic moments on the substrate in multiples of *dim_y*. Always use an even value (the code is not tested for odd values). The exact number of substrate moments in the y-direction can be calculated by (*substrate_dim_y* + 1) x *dim_y*. Note that the source code is only tested for systems where *substrat_dim_x* is equal to *substrat_dim_y*. The number of magnets on the substrate is not matched to the experiments but freely chosen.
- **angle_substrat:** Orientation angle of the substrate moments in the xz-plane. The value 0 corresponds to moments that align with the x-axis.
- **magnetic_moment_sub:** Magnitude of the magnetic moments on the substrate. For [1], this value is approximated by using the magnetization of the magnetic material used in the experiments and the volume of the cylinder-shaped substrate magnets.
- **tilt_substrat:** In theory, the slider and the substrate can be rotated relative to each other. This parameter is, however, not used in [1] and should be ignored by setting it to zero.
- **shift_x:** Initial displacement of the slider in the x-direction. The initial position of the slider is perfectly centered on the substrate (see the supporting information of [1]). The probe can be displaced from this position using *shift_x*.
- **shift_y:** Initial displacement of the slider in the y-direction. The initial position of the slider is perfectly centered on the substrate (see the supporting information of [1]). The probe can be displaced from this position using *shift_y*.
- **shift_z:** Separation between the substrate and the slider. This parameter gets overwritten by *scan_value*.
- **tilt_angle:** In theory, the slider and the substrate can be rotated relative to each other. This parameter is, however, not used in [1] and should be ignored by setting it to zero.
- **time_step:** Length of the simulation step.
- **max_time:** Maximum time of the simulations in seconds (time of the simulation, not real time). The simulation ends when it reaches *max_time*.
- **start_t_time:** Equilibration time. Until *start_t_time*, the slider idles at the initial position. After *start_t_time*, the simulation starts.
- **scan_value:** Interval of different slider-substrate separations for which the simulations should be performed.

[1] *Nonmonotonic Magnetic Friction from Collective Rotor Dynamics* by Hongri Gu, Anton Lüders, and Clemens Bechinger.
