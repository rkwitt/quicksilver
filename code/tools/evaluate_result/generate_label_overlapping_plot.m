close all, clear all;

% BEGINING USER INPUT
% input the directories + file names of the label overlapping result .mat 
% files generated from calculate_XXXX_overlap.m in this folder
LPBA_result_mat_dir = '';
IBSR_result_mat_dir = '';
CUMC_result_mat_dir = '';
MGH_result_mat_dir = '';
% END USER INPUT

%evaluate LPBA40 result
load ./Quicksilver_results/LPBA_results.mat;
load(LPBA_result_mat_dir);
results = [results, result_mean];
direc_name{end+1} = 'Your result';
figure(1), boxplot(results, 'labels', direc_name, 'PlotStyle', 'compact'), title('LPBA40')

%evaluate CUMC12 result
load ./Quicksilver_results/CUMC_results.mat;
load(CUMC_result_mat_dir);
results = [results, result_mean];
direc_name{end+1} = 'Your result';
figure(2), boxplot(results, 'labels', direc_name, 'PlotStyle', 'compact'), title('CUMC12')

%evaluate IBSR18 result
load ./Quicksilver_results/IBSR_results.mat;
load(IBSR_result_mat_dir);
results = [results, result_mean];
direc_name{end+1} = 'Your result';
figure(3), boxplot(results, 'labels', direc_name, 'PlotStyle', 'compact'), title('IBSR18')

%evaluate MGH10 result
load ./Quicksilver_results/MGH_results.mat;
load(MGH_result_mat_dir);
results = [results, result_mean];
direc_name{end+1} = 'Your result';
figure(4), boxplot(results, 'labels', direc_name, 'PlotStyle', 'compact'), title('MGH10')