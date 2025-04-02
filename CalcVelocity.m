function [velocity] = CalcVelocity(magnet_array_0, magnet_array_1, ...
    time_step)

    velocity = (magnet_array_1 - magnet_array_0) / time_step;
    
end