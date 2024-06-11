% Clear workspace and command window
clc;
clear;

% Define parameters
azimuth_angles = -90:.5:90; % Azimuth angles
elevation_angles = -90:.5:90;   % Elevation angles
[azimuth_grid, elevation_grid] = meshgrid(azimuth_angles, elevation_angles); % Create a grid of angles

sine_azimuth = sind(azimuth_grid); % Sine of azimuth angles
sine_elevation = sind(elevation_grid);     % Sine of elevation angles

speed_of_light = 300 / 10; % Speed of light
antenna_spacing = 15.4; % Antenna spacing
wave_number = 2 * pi / speed_of_light;
array_factor = zeros(size(azimuth_grid)); % Initialize Array Factor
scan_azimuth = 0; % Scan angle in azimuth
scan_elevation = 0;   % Scan angle in elevation
truncation_level = -40;
num_elements = 16; % Number of elements along one dimension (assuming square array)

% Calculate amplitude weights using Taylor taper for 1D
amplitude_1D = taylorTapperEven(num_elements / 2, 30, 5);

% Generate 2D amplitude weights for the square array
amplitude_2D = amplitude_1D' * amplitude_1D; % Outer product to form a 2D amplitude weight matrix

% Initialize phase array (assuming no phase taper)
phase_array = zeros(num_elements, num_elements);

% Calculate 2D Array Factor (AF) using nested loops
for m = 1:num_elements
    for n = 1:num_elements
        array_factor = array_factor + amplitude_2D(m, n) * exp(-1j * phase_array(m, n)) .* ...
             exp(-1j * wave_number * antenna_spacing * ((m - (num_elements+1)/2) * (sine_azimuth - sind(scan_azimuth)) + (n - (num_elements+1)/2) * (sine_elevation - sind(scan_elevation))));
    end
end

% Normalize and convert to dB
normalized_AF = abs(array_factor) ./ max(max(abs(array_factor)));
AF_dB = 20 * log10(abs(normalized_AF));
indices_below_trunc = find(AF_dB < truncation_level);
AF_dB(indices_below_trunc) = truncation_level;

% Plot the 2D Array Factor with colormap
figure;
colormap('jet'); % Use 'jet' colormap
mesh(azimuth_angles, elevation_angles, AF_dB);
xlabel('Azimuth Angle (degrees)');
ylabel('Elevation Angle (degrees)');
zlabel('Array Factor (dB)');
title('2D Antenna Array Factor');
grid on;
grid minor; % Add minor grid lines for better readability
colorbar; % Add color bar for reference
axis equal; % Set aspect ratio to equal

% Taylor Taper Function for Even Number of Elements
function [amplitude] = taylorTapperEven(N, SL, n)
    % Taylor patterns for discrete arrays
    % Refer IEEE transation on Antennas and Propagation
    % vol. AP-32, no. 10, October 1984. pg 1089-1093
    % by Alfred T. Villeneuve

    eta = 10^(SL / 20);
    Ele = 2 * N;
    mu0 = cosh((1 / (Ele - 1)) * log((eta + sqrt(eta * eta - 1))));
    sigma = (n * pi) / ((Ele) * acos((1 / mu0) * cos((2 * n - 1) * (pi / (2 * (Ele - 1))))));

    for cnt = 1:(Ele - 1)
        sai_p(cnt) = 2 * acos((1 / mu0) * cos((2 * cnt - 1) * (pi / (2 * (Ele - 1)))));
    end
    d_sai_p = sigma * sai_p;

    mul = 1;
    for q = 1:(n - 1)
        mul = mul * (sin(d_sai_p(q) / 2)) * (sin(d_sai_p(q) / 2));
    end
    numE0 = (Ele) * mul;

    mul = 1;
    for q = 1:(n - 1)
        mul = mul * sin((q * pi) / (Ele)) * sin((q * pi) / (Ele));
    end
    denE0 = mul;
    E0 = numE0 / denE0;

    for m = 1:(n - 1)
        num_mul = 1;
        for q = 1:(n - 1)
            num_mul = num_mul * sin(((m * pi) / (Ele)) - (d_sai_p(q) / 2)) * sin(((m * pi) / (Ele)) + (d_sai_p(q) / 2));
        end

        den_mul = 1;
        for q = 1:(n - 1)
            if q ~= m
                den_mul = den_mul * sin(((m - q) * pi) / (Ele)) * sin(((m + q) * pi) / (Ele));
            end
        end
        E(m) = ((Ele) * (-1)^m * num_mul) / (sin((m * pi) / (Ele)) * sin((2 * m * pi) / (Ele)) * den_mul);
    end

    for p = 1:N
        sum = 0;
        for m = 1:(n - 1)
            sum = sum + E(m) * cos(((p - 0.5) * m * pi) / (N));
        end
        Amp(p) = (1 / (Ele)) * (E0 + 2 * sum);
    end
    amplitude = [fliplr(Amp) Amp];
    amplitude = amplitude / max(amplitude);
end
