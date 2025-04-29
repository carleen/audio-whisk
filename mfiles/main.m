clear; clc; close all;

%% User-defined variables
% Get the files of interests
songloc = '/Users/boyercm1/Documents/misc/tempo_processing/audio_files/sparhawk_white-roses/*.wav';
ii = 2; % Index of song to use
show_plots = 1;

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
t_a = 0;
% Amount of time to look at
time_window = 2; % in seconds
% Number of downbeats to consider
N_D = 50;

end_time_seconds = round(max(t_bins));

bpm_ests_arr = zeros(100,1);
show_plots = 0;
frames = (1:time_window:end_time_seconds-time_window*2);

max_times = zeros(15,length(frames));
max_bpms = zeros(15, length(frames));
max_cc_vals = zeros(15, length(frames));

for ii=1:length(frames)
    t_a = frames(ii); % Start time
    [beat_locations] = tempoDownbeatLocation(sum_flux, t_bins, t_a, time_window, N_D, show_plots);
    max_times(:,ii) = beat_locations.cc_max_time_loc;
    max_bpms(:,ii) = beat_locations.best_bpm_from_normalized';
    max_cc_vals(:,ii) = beat_locations.cross_correlation_values;
end
