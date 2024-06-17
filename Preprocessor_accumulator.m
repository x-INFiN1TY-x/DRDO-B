clc; clear;

fileName_phase = 'phase_data.xlsx';
fileName_array_factor = 'array_factor_data.xlsx';
fileName_theta = 'theta_data.xlsx';

sheet_phase = 'Sheet1';
sheet_array_factor = 'Sheet1';
sheet_theta = 'Sheet1';

theta_data = readtable(fileName_theta, 'Sheet', sheet_theta);
theta_values = table2array(theta_data(1, :));

phase_data = readtable(fileName_phase, 'Sheet', sheet_phase);
phase_values = table2array(phase_data);

array_factor_data = readtable(fileName_array_factor, 'Sheet', sheet_array_factor);
array_factor_values = table2array(array_factor_data);

num_samples = size(phase_values, 1);

assert(size(phase_values, 2) == 18, 'Expected 18 columns for phase values');
assert(size(array_factor_values, 2) == 181, 'Expected 181 columns for array factor values');
assert(size(array_factor_values, 1) == num_samples, 'Expected the same number of rows for array factor values and phase values');
assert(length(theta_values) == 181, 'Expected 181 theta values');

beamwidth_deg = zeros(num_samples, 1);
peak_value_dB = zeros(num_samples, 1);
average_sll_dB = zeros(num_samples, 1);
rms_sll_dB = zeros(num_samples, 1);

for i = 1
array_factor_dB = array_factor_values(i, :);

scss
Copy code
normalized_array_factor = (abs(array_factor_dB) ./ max(abs(array_factor_dB)));
array_factor_dB = 20 * log10(abs(normalized_array_factor));

main_lobe_indices = find(array_factor_dB >= -3);
if ~isempty(main_lobe_indices)
    beamwidth_deg(i) = theta_values(max(main_lobe_indices)) - theta_values(min(main_lobe_indices));
else
    beamwidth_deg(i) = NaN;
end

peaks = findpeaks(normalized_array_factor);
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
end

results_table = table(phase_values, beamwidth_deg, peak_value_dB, average_sll_dB, rms_sll_dB);

writetable(results_table, 'calculated_metrics.xlsx');

save('phase_values.mat', 'phase_values');
save('beamwidth_deg.mat', 'beamwidth_deg');
save('peak_value_dB.mat', 'peak_value_dB');
save('average_sll_dB.mat', 'average_sll_dB');
save('rms_sll_dB.mat', 'rms_sll_dB');
