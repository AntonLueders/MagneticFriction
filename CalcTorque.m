% Calculates the total torque affecting the moments of the slider.
% This torque is needed for the Verlet algorithm to update the
% orientations.
function [torques] = CalcTorque(Angles, Velocities, ...
    friction_coefficient, pos_x, pos_y, substrat_array, substrat_x, ...
    substrat_y, shift_z, tilt_angle, magnetic_moment, ...
    magnetic_moment_sub)

    friction_torques = - friction_coefficient * Velocities;

    [n,m] = size(Angles);
    [k,l] = size(substrat_array);

    magnetic_torques = zeros(n,m);

    B_pre = 1.25663706*10^(-6) / (4 * pi);

    % Magnetic moments
    m_x = magnetic_moment * cos(Angles) * cos(tilt_angle);
    m_y = magnetic_moment * cos(Angles) * sin(tilt_angle);
    m_z = magnetic_moment * (-sin(Angles));

    m_x_substrat = magnetic_moment_sub * cos(substrat_array);
    m_y_substrat = zeros(k,l);
    m_z_substrat = magnetic_moment_sub * (- sin(substrat_array));
    
    for i = 1:n
        for j = 1:m

            % Magnetic field resulting from probe spins
            dis_x = pos_x(i,j)*ones(n,m) - pos_x;
            dis_y = pos_y(i,j)*ones(n,m) - pos_y;
            dis = (dis_x.^2 + dis_y.^2).^(1/2);
            dis(i,j) = 1;            

            m_dot_r = dis_x .* m_x + dis_y .* m_y;

            B_x = B_pre ...
                * ((3 * dis_x .* m_dot_r) ./ dis.^5 - m_x ./ dis.^3);
            B_y = B_pre ...
                * ((3 * dis_y .* m_dot_r) ./ dis.^5 - m_y ./ dis.^3);
            B_z = B_pre ...
                * (- m_z ./ dis.^3);

            B_x(i,j) = 0;
            B_y(i,j) = 0;
            B_z(i,j) = 0;

            % Magnetic fields resulting from substrat spins            
            dis_x_substrat = pos_x(i,j) * ones(k,l) - substrat_x;
            dis_y_substrat = pos_y(i,j) * ones(k,l) - substrat_y;
            dis_z_substrat = shift_z * ones(k,l);
            dis_substrat = (dis_x_substrat.^2  + dis_y_substrat.^2 + ...
                dis_z_substrat.^2).^(1/2);
            
            m_dot_r_substrat = dis_x_substrat .* m_x_substrat ...
                + dis_z_substrat .* m_z_substrat;

            B_x_substrat = B_pre ...
                * ((3 * dis_x_substrat .* m_dot_r_substrat) ...
                ./ dis_substrat.^5 - m_x_substrat ./ dis_substrat.^3);
            B_y_substrat = B_pre ...
                * ((3 * dis_y_substrat .* m_dot_r_substrat) ...
                ./ dis_substrat.^5 - m_y_substrat ./ dis_substrat.^3);
            B_z_substrat = B_pre ...
                * ((3 * dis_z_substrat .* m_dot_r_substrat) ...
                ./ dis_substrat.^5 - m_z_substrat ./ dis_substrat.^3);

            % Calculation of the total torque            
            torque_x = m_y(i,j) * (sum(B_z,'all') ...
                + sum(B_z_substrat,'all')) - m_z(i,j) * (sum(B_y,'all') ...
                + sum(B_y_substrat,'all'));
            torque_y = m_z(i,j) * (sum(B_x,'all') ...
                + sum(B_x_substrat,'all')) - m_x(i,j) * (sum(B_z,'all') ...
                + sum(B_z_substrat,'all'));

            magnetic_torques(i,j) = - sin(tilt_angle) * torque_x ...
                + cos(tilt_angle) * torque_y;
        end
    end
    
    torques = friction_torques + magnetic_torques;

end
