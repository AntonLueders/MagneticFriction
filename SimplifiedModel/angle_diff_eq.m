% Equations of motion of the simplified model
function dangledt = angle_diff_eq(t, angle, a, b, c, omega, gamma)

    % Substrate term
    % The system is "equilibrated" for t = 2*pi/omega
    if t < 2 * pi / omega
        coup1 = c * sin(angle(1)) / gamma;
        coup2 = c * sin(angle(2)) / gamma;
    else
        coup1 = c * sin(angle(1) - omega * (t - 2 * pi/omega)) / gamma;
        coup2 = c * sin(angle(2) - omega * (t - 2 * pi/omega)) / gamma;
    end
    
    % Self term + sublattice interaction term + substrate term
    dangledt = [- 2 * a * sin(angle(1)) * cos(angle(1)) / gamma ...
                + b * sin(angle(1) - angle(2)) / gamma + coup1; ...
                - 2 * a * sin(angle(2)) * cos(angle(2)) / gamma ...
                - b * sin(angle(1) - angle(2)) / gamma + coup2];
end
