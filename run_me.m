%% Description of run_me

% This script is segmented into several parts. First, the data (an
% example cell) is loaded. Then, 15 LN models are fit to the
% cell's spike train. Each model uses information about 
% position, head direction, running speed, theta phase,
% or some combination thereof, to predict a section of the
% spike train. Model fitting and model performance is computed through
% 10-fold cross-validation, and the minimization procedure is carried out
% through fminunc. Next, a forward-search procedure is
% implemented to find the simplest 'best' model describing this spike
% train. Following this, the firing rate tuning curves are computed, and
% these - along with the model-derived response profiles and the model
% performance and results of the selection procedure are plotted.

% Code as implemented in Hardcastle, Maheswaranthan, Ganguli, Giocomo,
% Neuron 2017
% V1: Kiah Hardcastle, March 16, 2017


%% Clear the workspace and load the data

clear all; close all; clc

% % load the data
% fprintf('(1/5) Loading data from example cell \n')
% load data_for_cell77

fprintf('(1/5) Loading data from cells \n')
load spikePos_openfield_D4.mat
load thetaEEG_JZ1_day4_nt11.mat
load JZ1posday4.mat
posStruct6 = pos{1,4}{1,6};
posStruct8 = pos{1,4}{1,8};
spikePosStruct = JZ1spikesPosD4;
thetaStruct6 = thetaEEG_JZ1_day4_epoch6_nt11;
thetaStruct8 = thetaEEG_JZ1_day4_epoch8_nt11;
numCells = length(spikePosStruct); % number of cells being analyzed 

% description of variables included:
% boxSize = length (in cm) of one side of the square box
% post = vector of time (seconds) at every 20 ms time bin
% spiketrain = vector of the # of spikes in each 20 ms time bin
% posx = x-position of left LED every 20 ms
% posx2 = x-position of right LED every 20 ms
% posx_c = x-position in middle of LEDs
% posy = y-position of left LED every 20 ms
% posy2 = y-posiiton of right LED every 20 ms
% posy_c = y-position in middle of LEDs
% filt_eeg = local field potential, filtered for theta frequency (4-12 Hz)
% eeg_sample_rate = sample rate of filt_eeg (250 Hz)
% sampleRate = sampling rate of neural data and behavioral variable (50Hz)

%% loop through all cells 
for p = 1:numCells
    clc
    
    %% fit the model
    fprintf('(2/5) Cell %d: Fitting all linear-nonlinear (LN) models\n',p)
    fit_all_ln_models

    %% find the simplest model that best describes the spike train
    fprintf('(3/5) Cell %d: Performing forward model selection\n',p)
    select_best_model
                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    %% Compute the firing-rate tuning curves
    fprintf('(4/5) Cell %d: Computing tuning curves\n',p)
    compute_all_tuning_curves

    %% plot the results
    fprintf('(5/5) Cell %d: Plotting performance and parameters\n',p) 
    plot_performance_and_parameters
    
end
