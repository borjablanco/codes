% fNIRS course 2022!
%
% Part II - Removal of systemic physiology in the preprocessed dc
%   Solution 1: Only Short-Channel Regression.
%   Soluiton 2: Only Additional Physiology.
%   Solution 3: Short-Channel Regression combined with Physiology. 
%   Solution 4: PCA: Regression of 1st Component.

% Clear environment
clear; close all; clc

% load preprocessed data from Part 1.
load Data_for_Part_II.mat;

% Solution 1: Short-Channel Regression
AdditionalRegressors = [];
[dc_only_SC,~] = PhysiologyRegression_GLM_fnirs_course...
(dc,SD,SSlist,AdditionalRegressors);

% Solution 2: Independent Physiological Measurements
AdditionalRegressors = Phys_data;
[dc_only_Phys,~] = PhysiologyRegression_GLM_fnirs_course...
(dc,SD,[],AdditionalRegressors);

% Solution 3: Short Channel Regression and Systemic Physiology
AdditionalRegressors = Phys_data;
[dc_SC_Phys,Stats] = PhysiologyRegression_GLM_fnirs_course...
(dc,SD,SSlist,AdditionalRegressors);

% Solution 4: Employ PCA to remove the 1st component
nSV = [1 1];
[dc_pca_one] = Perform_pca_regression_fnirs_course...
(dc,SD,nSV,BadChannels);

% Solution 5: Global signal regression
hbo = dc(:,:,1);
mhbo = mean(hbo,2);
hbr = dc(:,:,2);
mhbr = mean(hbr,2);

x_hbo = mhbo;
beta_data = pinv(x_hbo)*hbo;
reg_hbo = hbo - x_hbo*beta_data;

x_hbr = mhbr;
beta_data = pinv(x_hbr)*hbr;
reg_hbr = hbr - x_hbr*beta_data;
dc = cat(3, reg_hbo,reg_hbr, reg_hbo);

%%%%%% Plot correlation matrices of each case
plot_correlation_matrices_fnirs_course;
