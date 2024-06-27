clc;
clear;

fileName_phase = 'Phase.xlsx';
fileName_array_factor = 'AF_samples.xlsx';
fileName_theta = 'Theta.xlsx';

sheet_phase = 'Sheet1';
sheet_array_factor = 'Sheet1';
sheet_theta = 'Sheet1';

theta_data = readtable(fileName_theta, 'Sheet', sheet_theta);
theta_values = table2array(theta_data(1, :));

phase_data = readtable(fileName_phase, 'Sheet', sheet_phase);
phase_values = table2array(phase_data);

% Define array element coordinates (example)
N = size(phase_values, 2);
element_coordinates = rand(N, 2) * 10;  % Random coordinates within a 100x100 area

% Generate random (x, y) coordinates for evaluation points
num_samples = size(phase_values, 1);
x_coords = round(rand(num_samples, 1) * 100, 1);
y_coords = round(rand(num_samples, 1) * 100, 1);

% Initialize matrices for results
beamwidth_deg = zeros(num_samples, 1);
peak_value_dB = zeros(num_samples, 1);
average_sll_dB = zeros(num_samples, 1);
rms_sll_dB = zeros(num_samples, 1);
amplitudes_matrix = zeros(num_samples, N);  % Store amplitudes for each sample

for i = 1:num_samples
    x = x_coords(i);
    y = y_coords(i);
    
    % Generate random amplitudes between 0 and 1 with increments of 0.1
    amplitudes = randi([0, 10], 1, N) / 10;
    amplitudes_matrix(i, :) = amplitudes;
    
    % Calculate distances from array elements to point (x, y)
    r_n = sqrt((element_coordinates(:, 1) - x).^2 + (element_coordinates(:, 2) - y).^2);
    
    % Calculate array factor
    k = 2 * pi;  % Assuming wavelength = 1 for simplicity
    AF = sum(amplitudes .* exp(1j * (phase_values(i, :)' + k * r_n)));
    
    array_factor_dB = 20 * log10(abs(AF));
    
    % Ensure array_factor_dB has at least 3 elements
    if length(array_factor_dB) < 3
        disp('Insufficient data points for peak detection.');
        continue;
    end

    % Directly work with array factor without normalization
    array_factor_dB = 20 * log10(abs(AF));

    % Find main lobe indices and calculate beamwidth
    main_lobe_indices = find(array_factor_dB >= -3);
    if ~isempty(main_lobe_indices)
        beamwidth_deg(i) = theta_values(max(main_lobe_indices)) - theta_values(min(main_lobe_indices));
    else
        beamwidth_deg(i) = NaN;
    end

    % Check if array_factor_dB has at least 3 elements before calling findpeaks
    if length(array_factor_dB) >= 3
        peaks = findpeaks(array_factor_dB);
        sorted_peaks = sort(peaks, 'descend');
        if length(sorted_peaks) > 1
            peak_value_dB(i) = 20 * log10(sorted_peaks(2));
        else
            peak_value_dB(i) = 20 * log10(sorted_peaks(1));
        end

        if length(sorted_peaks) > 1
            sorted_peaks(1) = [];
            num_peaks = numel(sorted_peaks);
            sum_peaks = sum(sorted_peaks);
            average_sll_dB(i) = 20 * log10(sum_peaks / num_peaks);
        else
            average_sll_dB(i) = NaN;
        end

        if ~isempty(sorted_peaks)
            rms_sll = sqrt(mean(sorted_peaks .^ 2));
            rms_sll_dB(i) = 20 * log10(rms_sll);
        else
            rms_sll_dB(i) = NaN;
        end
    else
        peak_value_dB(i) = NaN;
        average_sll_dB(i) = NaN;
        rms_sll_dB(i) = NaN;
    end
end

results_table = table(x_coords, y_coords, amplitudes_matrix, phase_values, beamwidth_deg, peak_value_dB, average_sll_dB, rms_sll_dB);
writetable(results_table, 'calculated_metrics_with_coordinates.xlsx');

save('phase_values.mat', 'phase_values');
save('beamwidth_deg.mat', 'beamwidth_deg');
save('peak_value_dB.mat', 'peak_value_dB');
save('average_sll_dB.mat', 'average_sll_dB');
save('rms_sll_dB.mat', 'rms_sll_dB');
save('amplitudes_matrix.mat', 'amplitudes_matrix');
