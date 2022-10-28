function eq = aircraftPos(x)

    f = 1e8; % signal frequency (Hz)
    c = 3e8; % speed of light (m/s)
    lambda = c/f; % wavelength (m)
    dt = 0.1; % sampling period (s)

    xdot = 360 * (1000/3600); % aircraft speed (m/s)
    dx = xdot * dt; % change in X position (m)
    phi = [-33.1679 -33.1711 -33.1743]; % doppler shifts (Hz)
    rdot = zeros(1, length(phi)); % pseudorange rate preallocation

    for i = 1:length(phi)
    
    rdot(i) = -phi(i) * lambda; % pseudorange rate calculations

    end

    eq(1) = (( (x(1) + dx) * xdot ) / ( sqrt( (x(1) + dx)^2 + x(2)^2 ) )) - rdot(2);
    eq(2) = (( (x(1) + 2*dx) * xdot ) / ( sqrt( (x(1) + 2*dx)^2 + x(2)^2 ) )) - rdot(3);
    
end