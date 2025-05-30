% This function updates the spatial position of the slider relative to
% the substrate.
function [pos_in_grid_x_temp, pos_in_grid_y_temp] ...
            = Translate(pos_in_grid_x, pos_in_grid_y, v_x, v_y, ...
            time_step)

    pos_in_grid_x_temp = pos_in_grid_x + v_x * time_step;
    pos_in_grid_y_temp = pos_in_grid_y + v_y * time_step;

end
