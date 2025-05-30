% Determines the substrate torque acting on each of the moments of
% the slider and saves it as a matrix
function T = CalcTorqueGrid(Angles, substrat_array, tilt_angle, ...
    magnetic_moment, magnetic_moment_sub, shift_z, pos_x, pos_y, ...
    substrat_x, substrat_y)
    
    [n,m] = size(Angles);
    [k,l] = size(substrat_array);

    Bx = 0;
    By = 0;
    Bz = 0;

    T = zeros(n,m);

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

            Bx = Bx + sum(B_x_substrat,'all');
            By = By + sum(B_y_substrat,'all');
            Bz = Bz + sum(B_z_substrat,'all');

            torque_x = m_y(i,j) * sum(B_z_substrat,'all')...
               - m_z(i,j) * sum(B_y_substrat,'all');
            torque_y = m_z(i,j) * sum(B_x_substrat,'all')...
               - m_x(i,j) * sum(B_z_substrat,'all');

            T(i,j) = - sin(tilt_angle) * torque_x + cos(tilt_angle) * torque_y;
        end
    end    
end
