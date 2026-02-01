% Calculates the magnetic interactions between the slider and the substrate
% This means that the magnetic friction is calculated in this function
function [Fx, Fy, Fz] = CalcForce(Angles, substrat_array, tilt_angle, ...
    magnetic_moment, magnetic_moment_sub, shift_z, pos_x, pos_y, ...
    substrat_x, substrat_y)
    
    [n,m] = size(Angles);
    [k,l] = size(substrat_array);

    magnetic_forces_x = zeros(n,m);
    magnetic_forces_y = zeros(n,m);
    magnetic_forces_z = zeros(n,m);

    F_pre = 1.25663706*10^(-6) * 3 / (4 * pi);

    % Magnetic moments
    m_x = magnetic_moment * cos(Angles) * cos(tilt_angle);
    m_y = magnetic_moment * cos(Angles) * sin(tilt_angle);
    m_z = magnetic_moment * (-sin(Angles));

    m_x_substrat = magnetic_moment_sub * cos(substrat_array);
    m_y_substrat = zeros(k,l);
    m_z_substrat = magnetic_moment_sub * (- sin(substrat_array));

    for i = 1:n
        for j = 1:m

            dis_x = pos_x(i,j)*ones(n,m) - pos_x;
            dis_y = pos_y(i,j)*ones(n,m) - pos_y;
            dis_z = zeros(n,m);
            dis = (dis_x.^2 + dis_y.^2).^(1/2);
            dis(i,j) = 1;

            m1_dot_r = m_x .* dis_x + m_y .* dis_y;

            m2_dot_r = m_x(i,j) * dis_x + m_y(i,j) * dis_y;

            m1_dot_m2 = m_x(i,j) * m_x + m_y(i,j) * m_y;

            % Does not contribute due to Newtons 3. law. 
            % Adds to zero.
            F_x = F_pre ./ dis.^5 .* (m1_dot_r .* m_x(i,j) ...
                + m2_dot_r .* m_x ...
                + m1_dot_m2 .* dis_x ....
                - 5 * m1_dot_r .* m2_dot_r ./ dis.^2 ...
                .* dis_x);
            F_y = F_pre ./ dis.^5 .* (m1_dot_r .* m_y(i,j) ...
                + m2_dot_r .* m_y ...
                + m1_dot_m2 .* dis_y ....
                - 5 * m1_dot_r .* m2_dot_r ./ dis.^2 ...
                .* dis_y);
            F_z = F_pre ./ dis.^5 .* (m1_dot_r .* m_z(i,j) ...
                + m2_dot_r .* m_z ...
                + m1_dot_m2 .* dis_z ....
                - 5 * m1_dot_r .* m2_dot_r ./ dis.^2 ...
                .* dis_z);

            F_x(i,j) = 0;
            F_y(i,j) = 0;
            F_z(i,j) = 0;

            % Calculation of distance vector            
            dis_x_substrat = pos_x(i,j) * ones(k,l) - substrat_x;
            dis_y_substrat = pos_y(i,j) * ones(k,l) - substrat_y;
            dis_z_substrat = shift_z * ones(k,l);
            dis_substrat = (dis_x_substrat.^2  + dis_y_substrat.^2 + ...
                dis_z_substrat.^2).^(1/2);

            % Calculation of the dot products in the force equation
            m1_dot_r = m_x_substrat .* dis_x_substrat ...
                + m_y_substrat .* dis_y_substrat ...
                + m_z_substrat .* dis_z_substrat;

            m2_dot_r = m_x(i,j) * dis_x_substrat ...
                + m_y(i,j) * dis_y_substrat ...
                + m_z(i,j) * dis_z_substrat;

            m1_dot_m2 = m_x(i,j) * m_x_substrat ...
                + m_y(i,j) * m_y_substrat ...
                + m_z(i,j) * m_z_substrat;

            % Each force acting on spin i,j
            Fx_all = F_pre ./ dis_substrat.^5 .* (m1_dot_r .* m_x(i,j) ...
                + m2_dot_r .* m_x_substrat ...
                + m1_dot_m2 .* dis_x_substrat ....
                - 5 * m1_dot_r .* m2_dot_r ./ dis_substrat.^2 ...
                .* dis_x_substrat);
            Fy_all = F_pre ./ dis_substrat.^5 .* (m1_dot_r .* m_y(i,j) ...
                + m2_dot_r .* m_y_substrat ...
                + m1_dot_m2 .* dis_y_substrat ....
                - 5 * m1_dot_r .* m2_dot_r ./ dis_substrat.^2 ...
                .* dis_y_substrat);
            Fz_all = F_pre ./ dis_substrat.^5 .* (m1_dot_r .* m_z(i,j) ...
                + m2_dot_r .* m_z_substrat ...
                + m1_dot_m2 .* dis_z_substrat ....
                - 5 * m1_dot_r .* m2_dot_r ./ dis_substrat.^2 ...
                .* dis_z_substrat);

            % Matrix that contains the total forces acting on the spins
            magnetic_forces_x(i,j) = sum(Fx_all,'all') + sum(F_x,'all');
            magnetic_forces_y(i,j) = sum(Fy_all,'all') + sum(F_y,'all');
            magnetic_forces_z(i,j) = sum(Fz_all,'all') + sum(F_z,'all');

        end
    end

    % Total force acting on the probe
    Fx = sum(magnetic_forces_x, 'all');
    Fy = sum(magnetic_forces_y, 'all');
    Fz = sum(magnetic_forces_z, 'all');
    
end

