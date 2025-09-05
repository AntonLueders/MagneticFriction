# Introduction

MagneticFriction is a many-body simulation software that calculates the magnetic friction resulting from the relative motion of two periodic magnetic layers (one with moments that can rotate around an axis perpendicular to the movement direction). It utilizes the Verlet algorithm known from Molecular dynamics simulations in combination with a Stokes law model for the single particle dissipation and integrates the rotations of the moments of the slider, which is the layer with the moveable moments. MagneticFriction is a part of the supporting information of 

[1] *Nonmonotonic Magnetic Friction from Collective Rotor Dynamics* by Hongri Gu, Anton Lüders, and Clemens Bechinger.

The simulation software reproduces the qualitative and quantitative results of the experiments of [1] remarkably well.  All information and methods employed by MagneticFriction are described in detail in the supporting information of [1]. While MagneticFriction can be used to repeat the corresponding simulations, we strongly encourage a full independent reproduction using the information in said supplemental information. 

# Table of content:
 - [Requirements](#Requirements)
 - [How to use MagneticFriction](#Use)
 - [What MagneticFriction does](#What)
 - [How to interprete the generated data](#Data)
 - [Assumptions and simplifications](#Assumptions)
 - [Contents of the particular files](#Contents)
 - [Further information and examples](#Info)
 - [Disclaimer](#Disclaimer)
 - [Additional folders independent of MagneticFriction](#Simplified)
 - [Cite this software](#Cite)
 
 <a id="Requirements"></a>
# Requirements

The simulation software is written as a Matlab script, and it can be started as usual from within Matlab. The code was tested with **MatlabR2024b**, which was also used to generate the data given in [1].
 
 <a id="Use"></a>
# How to use MagneticFriction

To run the simulation software, all files must be in the same folder. The main script is **MagneticArray.m**, which can be started via Matlab. The parameters for the simulation can be adjusted at the top of **MagneticArray.m**. The ones used to generate the data for [1] are given below as an example.

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
scan_value=0.0050:0.0001:0.013; % Varied paramter of the simulations

for scan=1:length(scan_value)
    
    shift_z=scan_value(scan);
```
.
.
.

**Notes:** The parameters in the example are matched to the experiments of [1] (the shaft friction and the number of moments on the substrate are chosen freely). Also, *substrat_dim_x* and *substrate_dim_y* are given in units of the slider array size. They control the number of additional layers of arrays of moments that are added to the substrate. A value of 4 means that the substrate has a width of 5 slider-sized layers of moments, where the slider is centered in the middle at the start of the simulation. Please only use even numbers for these parameters.
 
 <a id="What"></a>
 # What MagneticFriction does

Detailed information about the simulation software is given in the Methods section and the supporting information of [1]. Using the defined parameters, the script generates a rigid array of magnetic moments (i.e., the slider) that hovers over a substrate of magnetic moments. Here, the magnetic moments of the slider can rotate in the xz-plane while the moments of the substrate are fixed. The dynamics of the moments of the slider are updated using a Molecular dynamics simulation scheme (Verlet algorithm) with an additional dissipation term, which models the shaft friction. 

After an equilibration time (defined by the user), the slider starts to translate over the substrate. Because of the magnetic interactions between the moments of the slider and between the slider and the substrate, magnetic torques act on the moments, and they start to rotate (in complex patterns) while translating. The magnetic interactions between the moments are modeled via standard dipole interactions. During the sliding, the simulation code calculates observables, such as the total interaction force between the slider and the substrate. The total interaction force can be averaged to obtain the finite magnetic friction. The calculations are repeated for multiple different separations between the slider and the substrate.

<a id="Data"></a> 
 # How to interprete the generated data

All data are given in SI units. For each separation between the slider and the substrate, the simulation software generates a **.mat** file, which stores the corresponding data. Given a predefined sampling rate (**Force_calc_rate** in **MagneticArray.m**), the codes stores for each sampling frame 

- the corresponding simulation time (**time**),
- the current orientation of the moments (**array**),
- the magnetic interaction force between the slider and the substrate (**force** - x-component, y-component, z-component),
- the particular subsrate torques acting on the different moments of the slider (**torque_grid**),
- all information needed to extract the sum of the torques and the torque hysteresis (**hysteresis** - sum of torques acting on moments of sublattice one, sum of torque acting on sublattice two, x-component average magnetic field, y-component average magnetic field, z-component average magnetic field, x-component of the magnetization of sublattice one, y-component of the magnetization of sublattice one, z-component of the magnetization of sublattice one, x-component of the magnetization of sublattice two, y-component of the magnetization of sublattice two, z-component of the magnetization of sublattice two; see [1] for information on how the slider is divided in two sublattices and the simplified model),
- and the energy terms of the system (**energy** - internal magnetic energy, total interaction energy between the slider and the substrate).
 
  <a id="Assumptions"></a>
# Assumptions and simplifications

Compared to the experimental setup of [1], MagneticFriction utilizes the following approximations and simplifications:

- **Dipole interactions** All magnetic interactions between the cylinder-shaped magnets are modeled via classical dipole interactions.
- **Stokes-like shaft friction** For the dissipation term $` f `$ that enters the equations of motion of the particular moments of the slider, a Stokes law $` f(\omega) = - \gamma \omega `$ is utilized.
- **Angular velocity** The angular velocities of the particular moments are approximated via a simple difference quotient. Note that this angular velocity enters the integration step through the Stokes friction law that models the shaft friction. 

For a detailed description of the limitations of the simulations, see the supporting information of [1]. Note that MagneticFriction generates data with excellent qualitative and quantitative agreement with the experiments of [1].

 
 <a id="Contents"></a>
# Contents of the particular files

MagneticFriction consists of multiple Matlab functions stored in their individual files. Here, a summary of the corresponding content of the particular files is given.

- **MagneticArray.m:** Main file which contains the simulation loop. It calls functions stored in the other files. Here, the parameters are adjusted, and a simulation is started by running this script.
- **CalcVelocity.m** Approximates the angular velocities of the moments, which are later used to calculate the shaft friction. 
- **CalcTorque.m** Calculates the total torque (shaft friction, internal magnetic interactions of the slider, and torques resulting from the substrate). These torques are used to update the angles of the moments in the integration step.
- **Move.m** Integration step based on the Verlet algorithm.
- **Translate.m** Moves the slider relative to the substrate.
- **CalcForce.m** Calculates the total magnetic forces acting on the slider. The average of the x-component of this force is the magnetic friction.
- **CalcEnegery.m** Determines the total internal interaction energy of the slider and the total energy of the interactions between the slider and the substrate.
- **CalcHysteresis.m** Calculates the sum of the torques affecting the moments due to the substrate, the magnetization of the slider, and the average substrate field.
- **CalcTorqueGrid.m** Evaluates and stores the total torque acting on each moment of the slider separately. 

To understand how MagneticFriction works, it is suggested to start by reading the main simulation loop in **MagneticArray.m**. If details for the specific function used in the main loop are needed, it is then possible to navigate to the correct subfile that contains the particular function.

 <a id="Info"></a>
# Further information and examples

A detailed explanation of the parameters needed to run the simulations are given in [Parameters](docs/Parameter.md).

Pseudocode for MagneticFriction can be found in [Pseudocode](docs/Pseudocode.md).

Please read the **Method section** and **Section S14** of the supporting information of [1] for information on the employed model, the simulation algorithms, and the physics behind **MagenticFriction**. 

 <a id="Disclaimer"></a>
# Disclaimer

Note that the presented code is **not** an ongoing software project with active support. Instead, it is part of the supplemental information of [1] and reflects the state of the software used to obtain the corresponding numerical results. For any questions regarding [1] and the corresponding data, please contact the corresponding authors.

 <a id="Simplified"></a>
# Additional folders independent of MagneticFriction

This repository also contains the folder **SimplifiedModel** with Matlab code to solve the simplified model (as introduced in [1]) numerically. These files are not related to the main simulation code. For questions regarding the numerical solution of the simplified model, please contact the corresponding authors.

 <a id="Cite"></a>
# Cite this software

If you use this simulation software or parts of it, please cite

[1] *Nonmonotonic Magnetic Friction from Collective Rotor Dynamics* by Hongri Gu, Anton Lüders, and Clemens Bechinger.
