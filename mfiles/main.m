clear; clc; close all;

%% User-defined variables
% Get the files of interests
songloc = '/Users/boyercm1/Documents/misc/tempo_processing/audio_files/sparhawk_white-roses/*.wav';
ii = 2; % Index of song to use
show_plots = 0;

nfft = 512; % Gives time resolution of ~ 10ms
noverlap = 0; % No overlap needed
min_thresh = -60;

f_min = 100; % Hz
f_max = 10e3; % Hz



%% File I/O
songs = dir(songloc);
filename = fullfile(songs(ii).folder, songs(ii).name);


%% Frequency-domain processing
[spectrogram_struct] = songSpectrogramProcess(filename, nfft, noverlap, ...
    min_thresh, show_plots);

%% Energy Flux Calculation
% Calculate the Energy Flux (Laroche, J. "Efficent Tempo and Beat Tracking
% in Audio Recordings")
show_plots = 1;
t_bins = spectrogram_struct(1).t;
f_bins = spectrogram_struct(1).f;
[spectrogram_struct, sum_flux] = fluxCalculator(spectrogram_struct, ...
    t_bins, f_bins, f_min, f_max, show_plots);

%% Tempo and Downbeat Location Estimation
% Start Time
t_a = 50;
% Amount of time to look at
time_window = 2;
% Number of downbeats to consider
N_D = 100;

bpm_ests_arr = zeros(100,1);
show_plots = 0;
for ii=1:100
    t_a = ii; % Start time
    [beat_locations] = tempoDownbeatLocation(sum_flux, t_bins, t_a, time_window, N_D, show_plots);
    bpm_est = beat_locations.bpm_candidates;
    upper = bpm_est(bpm_est>100);
    bpm_ests_arr(ii) = upper(1);
end


%% TO DO:
% For each candidate tempos, the 10-15 best candidate downbeat locations
% are selected






