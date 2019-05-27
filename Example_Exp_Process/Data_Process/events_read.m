%% Paths and Filenames
clear; clc; close all; 
dpath = '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/cat_004/meas/abr_powline/abr_run';
numrun = 5; raw = 2;

for i= 1:numrun
    %% Read Data
    if raw == 1
      fname = strcat(dpath, num2str(i),'/ABR_run',num2str(i),'_raw.fif');
      ename = strcat(dpath, num2str(i),'/ABR_run',num2str(i),'_raw-eve.fif');
      %etname = strcat(dpath,num2str(i),'/ABR_run',num2str(i),'_raw.eve');
    else
      fname = strcat(dpath, num2str(i),'/ABR_run',num2str(i),'_raw-1.fif');
      ename = strcat(dpath, num2str(i),'/ABR_run',num2str(i),'_raw-1-eve.fif');
      %etname = strcat(dpath,num2str(i),'/ABR_run',num2str(i),'_raw-1.eve');
    end
    
    %% Direct Making of Event File
    include{1} = 'MISC006'; threshold = 0.05;
    delshift = round(9.775/1000*5000);
    mne_make_combined_event_file3(fname,ename,include,0,threshold,delshift,i);
    eventlist = mne_read_events(ename); %size(eventlist)
    
end