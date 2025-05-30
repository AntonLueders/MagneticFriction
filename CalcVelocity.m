% Approximates the angular velocity of the moments of the slider. The
% result is given to "CalcTorque" for the calculation of the shaft
% friction.
function [velocity] = CalcVelocity(magnet_array_0, magnet_array_1, ...
    time_step)

    velocity = (magnet_array_1 - magnet_array_0) / time_step;    
end
