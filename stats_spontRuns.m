function [ total_mean_rate, s_error, total_median_rate, IQR, CI_95] = stats_spontRuns( analysis_dir, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% ***IMPORTANT*** Requires creation of an "analysis directory" inside the cell directory of
% interest.  All rep directories to be included in the analysis must be
% moved to the analysis directory prior to running this function.

% Currently, you need to be in the analysis_dir folder to run this.  To be
% updated.
%
% VARARGIN can be used to define an array list containing the directory
% numbers to be analyzed; this allows one or a few of the subdirectories to be analyzed
% while excluding the remainder.  (One small problem with this code: the
% index "i" may number the subdirectories arbitrarily, i.e. not mirroring
% the numberings in the directory names.  NEED TO TEST THIS.)

reps_dir = dir(analysis_dir);
rate_array = [];
matlab_root = pwd;

if nargin < 2,
    to_analyze = (1:length(reps_dir));
else
    to_analyze = varargin;
end

for i=1:length(reps_dir),
    if any(ismember(to_analyze,i)),
        plotit = 0;
        [f_rate] = analyze_spontRate([matlab_root filesep reps_dir(i).name], plotit);
        rate_array(end+1,1) = f_rate;
        fprintf('%s \n',['Analyzed folder ', reps_dir(i).name]);
    end
end

total_mean_rate = mean(rate_array(:,1));
s_error = std2(rate_array(:,1))/(sqrt(length(rate_array(:,1))));
total_median_rate = median(rate_array(:,1));
low_IQR = total_median_rate - prctile(rate_array(:,1),25);
high_IQR = total_median_rate + prctile(rate_array(:,1),75);
IQR = [low_IQR;high_IQR];
IQR_length = high_IQR - low_IQR;
low_CI = total_median_rate - (1.58 * (IQR_length/sqrt(length(rate_array(:,1)))));
high_CI = total_median_rate + (1.58 * (IQR_length/sqrt(length(rate_array(:,1)))));
CI_95 = [low_CI;high_CI];
    
end

