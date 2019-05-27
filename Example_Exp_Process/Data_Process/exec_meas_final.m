%% Initialization and Parameters
clear; clc; close all;
prefix = 'cat_004';                                                         %Subject name
meastype = 'meg_eeg';                                                       %MEG or EEG forward solution
datatype = 'AEP';                                                           %Paradigm of data: Illustrations, SEP, AEP
datapath = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/', datatype, '/',prefix);
plot_on = 1;                                                                %Check the various covariances and whiteners
saveon = 0;                                                                 %Save the covariance and whitener information

%% Measurement Process: Select Channels, Derive Noise Covariances, Projectors, Whitener
append = '_MLR';                                                            %MLR: indicator for diary name, meas_fname
cat004_measinfo_megeeg(prefix, datatype, meastype, datapath, plot_on, append, saveon);%Derive whitener and noise covariances - for use in srcfwd and inv solns

append = '_HFABR';                                                          %HFABR: indicator for diary name, meas_fname
cat004_measinfo_megeeg(prefix, datatype, meastype, datapath, plot_on, append, saveon);%Derive whitener and noise covariances - for use in srcfwd and inv solns
