% fNIRS course 2022!
%
% Part I - Preprocessing required before the
% removal of systemic physiology
%   1 - Check the quality of the channels;
%   2 - Motion Corretion in the Optical Density;
%   3 - Compute hemoglobin concentration changes (dc);
%   4 - Band Pass filter dc and additional physiological measurements if
%   available.
%   5 - Remove border effects.

% Clear environment
clear

% Load data from one participant
load('Data_for_Part_I.mat')

% Get Light Intensity, SD, and additional physiological measurements
d = data.d;
SD = data.SD;

Phys_data = [data.Phys.MAP_d,...
    data.Phys.HR_d,...
    data.Phys.CapData];

% ***********************************

% List of short channels for the used probe
SSlist = [8 29 52 66 75 92 112 125];

% Find channels with low SNR
% The list of bad-channels will be used latter on further analysis
BadChannels = MarkBadChannels_fnirs_course(d,SD);

% Compute Optical Density
dOD = hmrIntensity2OD(d);

%*** Motion Correction
% Spline interpolation followed by wavelet decomposition
%dOD = Hybrid_motion_correction(dOD,SD);

% Compute Hemoglobin Concentration changes
dc = hmrOD2Conc...
    (dOD, SD, [6 6 6 6]);

% Permute dc
dc = permute(dc,[1 3 2]);

% Band-Pass Filter Hemoglobin concentrations
dc = hmrBandpassFilt...
    (dc, SD.f, 0.009, 0.08);

% Band-Pass Filter Additional Physiological Measurements
Phys_data = hmrBandpassFilt...
    (Phys_data, SD.f, 0.009, 0.08);

% Remove border effects
dc = dc(200:end-200,:,:);
Phys_data = Phys_data(200:end-200,:);

fc_hbo = corr(squeeze(dc(:,:,1)));

% Remove Autocorrelation 
[pw_dc, P_all] = RemoveAutocorrelation_dc_fnirs_course...
    (dc,SD);

% Plot autocorrelation function before and after PW
figure
data_autocorr = xcorr(squeeze(pw_dc(:,1,1)), 'coeff');
plot(data_autocorr(2339:end))
hold on
data_hbo = xcorr(squeeze(dc(:,1,1)), 'coeff');
plot(data_hbo(2417:end))
legend('pw', 'hbo')
fc_hbo_aut = corr(squeeze(pw_dc(:,:,1)));

% Plot connectivity matrices after PW
figure
subplot(131)
imagesc(fc_hbo, [-1 1])
subplot(132)
imagesc(fc_hbo_aut, [-1 1])
subplot(133)
imagesc(fc_hbo-fc_hbo_aut, [-1 1])
colormap jet

% Plot PSD of data
figure
subplot(131)
[pxx1,f] = pwelch(squeeze(dc(:,:,1)),[],[],[],SD.f);
plot(f,10*log10(pxx1))
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
xlim([0 1])
subplot(132)
[pxx2,f] = pwelch(squeeze(pw_dc(:,:,1)),[],[],[],SD.f);
plot(f,10*log10(pxx2))
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
xlim([0 1])
subplot(133)
plot(f,10*log10(mean(pxx1,2)))
hold on
plot(f,10*log10(mean(pxx2,2)))
xlim([0 1])
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
legend('hbo', 'hbo_pw', 'interpreter', 'none')

% Compute Pearson Correlation Coefficient
CorrMatrix = ...
    Compute_correlation_coefficient_fnirs_course...
    (pw_dc,BadChannels);

% Plot Sensory Motor Network
plot_seed_based_sphere_style_fnirs_course...
    (CorrMatrix(:,:,3),BadChannels,[-1 1]);


% QUESTIONS
% Why data is not white after prewhitening? figure 1
% Connectivity matrices very similar before and after pw. figure 2
% Why hbo and hbr same order?

% Save Data
%  save('Data_for_Part_II',...
%      'dc','SD','Phys_data','SSlist','BadChannels');





