close all, clear all
addpath(genpath('./medical_image_matlab_io'))
load g_Labels.mat
result = zeros(10*9, length(Labels));


%BEGIN USER DEFINED INPUT. CHANGE THESE DIRECTORIES BEFORE RUNNING THE CODE
% user defined directory of the label files (suppose are in the same directory) and registration result directory (one directory for one phiinv.mhd files)
label_files_dir = ''
registration_results_dir = cell(10, 10)
for from = 1:10
	for to = 1:10
		if (from == to)
			continue
		end
		registration_results_dir{i}{j} = ''
	end
end
% the directory and name of the output.mat file
output_name = ''
%END USER DEFINED INPUT

label_images = cell(10, 1);
for i = 1:10
	label_images{i} = load_nii(strcat(label_files_dir, 'g', num2str(i), '.nii'));
end

base_idx = 1;
for from = 1:10
	for to = 1:10
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

result_mean = zeros(10*9, 1);
for i = 1:10*9
	single_result = result(i, :);
	idx = find(~isnan(single_result));
	result_mean(i) = mean(single_result(idx));
end


save(output_name, 'result_mean');

		