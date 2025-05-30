% Integration step of the equations of motion of the slider spins. 
% The code is based on the Verlet algorithm known from MD simulations.
function [angle_0, angle_1] = Move(magnet_array_0, magnet_array_1, ...
    torques, time_step, moment_of_inertia)

    angle_1 = 2 * magnet_array_1 - magnet_array_0 ...
        + torques / moment_of_inertia * time_step^2;
    angle_0 = magnet_array_1;
    
end
