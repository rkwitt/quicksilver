close all, clear all
addpath(genpath('./medical_image_matlab_io'))
load c_Labels.mat
result = zeros(18*17, length(Labels));

%BEGIN USER DEFINED INPUT. CHANGE THESE DIRECTORIES BEFORE RUNNING THE CODE
% user defined directory of the label files (suppose are in the same directory) and registration result directory (one directory for one phiinv.mhd files)
label_files_dir = ''
registration_results_dir = cell(18, 18)
for from = 1:18
	for to = 1:18
		if (from == to)
			continue
		end
		registration_results_dir{i}{j} = ''
	end
end
% the directory and name of the output.mat file
output_name = ''
%END USER DEFINED INPUT

% load data
label_images = cell(18, 1);
for i = 1:18
	label_images{i} = load_nii(strcat(label_files_dir, 'c', num2str(i), '.nii'));
end

base_idx = 1;
for from = 1:18
	for to = 1:18
		if from == to 
			continue;
		end
		label_from = label_images{from}.img;
		label_to = label_images{to}.img;
		phi = loadMETA(strcat(registration_results_dir{from}{to}, 'phiinv.mhd'));
		warp_result = wrap_image_nn(label_from, phi);
		for label_idx = 1:length(Labels)
			warp_idx = find(warp_result(:) == Labels(label_idx));
			to_idx = find(label_to(:) == Labels(label_idx));
			result(base_idx, label_idx) = length(intersect(warp_idx, to_idx)) / length(to_idx);
		end
		single_result = result(base_idx, :);
		idx = find(~isnan(single_result));
		mean(single_result(idx))
		base_idx = base_idx+1;
	end
end

result_mean = zeros(18*17, 1);
for i = 1:18*17
	single_result = result(i, :);
	idx = find(~isnan(single_result));
	result_mean(i) = mean(single_result(idx));
end


save(output_name, 'result_mean');

		