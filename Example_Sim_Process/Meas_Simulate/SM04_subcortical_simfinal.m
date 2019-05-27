%% INITIALIZE
% General Paths and Parameters
clear; clc; close all;
prefix = 'SM04'; datatype = 'SEP'; meastype = 'meg';                        %Subject and Data Characteristics
add_allpaths;                                                               %Add Paths
overlap = '_overlap';    %overlap = '_overlap' or '_no-overlap';            %indicator re. overlap
diary(strcat(prefix,'_',meastype,'_sim_whitener_matlog',overlap));          %log of simulation data

% Simulation Parameters
saveon = 0; 
savename = strcat(datapath, 'meas/', prefix,'_vpl_somato_meg_sim',overlap); %Simulated Data Saving Parameters
sfreq = 3000;                                                               %sampling frequency [Hz]
scort_t_time = 15;                                                          %time [msec] for subcortical simulation
cort_t_time = 255;                                                          %time [msec] for cortical simulation 
tot_t_time = scort_t_time + cort_t_time;                                    %total time for all simulations [msec]
cort_repeat_delay = 0; %125; %floor(cort_t_time/num_cort_repeats);          %time [msec] b/w cortical timecourse repeats (0 = no delay)
if strcmp(overlap,'_overlap')
    num_scort_repeats = 10;                                                 %number of times to repeat subcortical time course - equalize norm to cort
    scort_repeat_delay = 25;                                                %time [msec] b/w subcortical timecourse repeats (0 = no delay)
else
    num_scort_repeats = 1;                                                  %subcortical time course only once - SUPP figure
end
SNR = 5;                                                                    %signal to noise ratio - 7 dB
ploton = 0;  savemem = 0;  p_sing = 6; nthresh = 95;                        %parameters for reduced SVD (for VPL)
t0 = 221;                                                                   %time point when response starts in meas evoked file 
trim_end = 18;                                                              %take out last trim_end points to avoid weird behaviors 
ico_num = 2;                                                                %cortical ico to consider

% Load Relevant Meas, Src and Fwd Files
meas_fname = strcat(datapath, 'meas/', prefix,'_SomSens_SingleNoTask_offl_SSP.fif'); %measurement file
meas = fiff_read_evoked_all(meas_fname);                                    %measurement file
src_fname = strcat(datapath, 'src/', prefix,'-ico-2p-src.fif');             %cortical source space - ico-2
src = mne_read_source_spaces(src_fname);                                    %cortical source space file - ico-2
fwd_fname = strcat(datapath, 'fwd/', prefix,'-meg-all-fixed-fwd.fif');      %fwd file - uses 3-BEM, no mindist
fwd = mne_read_forward_solution(fwd_fname);                                 %cortical forward file
fwd4_fname = strcat(datapath, 'fwd/', prefix,'-meg-ico-4-fixed-fwd.fif');   %fwd file - uses 3-BEM, no mindist
fwd4 = mne_read_forward_solution(fwd4_fname);                               %cortical forward file - ico-4
ssvds_fname = strcat(datapath, 'fwd/140620_Thesis_SrcSpace/',prefix,'_subc_meg_volumes.mat'); %subcortical srcpace (Jun2014 manual VPL)
load(ssvds_fname,'sdiv_fwd');                                               %subcortical divisions file with VPL for simulation
curr_den_fname = strcat(datapath,'fwd/','curr_den.mat');                    %file containing regional current densities
load(curr_den_fname,'cort_cd','th_cd');                                     %cortical and thalamic current densities

% Labels for Regions to Simulate
X(1).label_name = 'cS1';    X(1).label = '058397-lh.label';                 %ICO-2 label cS1
X(2).label_name = 'cS2';    X(2).label = '080245-lh.label';                 %ICO-2 label cS2
X(3).label_name = 'PPC';    X(3).label = '028238-lh.label';                 %ICO-2 label PPC: also 028238, 042471, 047227, 028022, 041027
X(4).label_name = 'iS2';    X(4).label = '084729-rh.label';                 %ICO-2 label iS2
label_string = {X(:).label};                                                %label listing by region
for i = 1:length(sdiv_fwd)
    all_regions{i} = sdiv_fwd(i).reg_name;
end
vpl_sdiv = find(~cellfun(@isempty,strfind(all_regions,'lt')),1);            %Division Corresponding to VPL

% Initialize Time Vector and STC Files
N = size(fwd.sol.data,1);                                                   %# sensors
for ii = 1:2
    stc(ii).tmin = 0;                                                       %initialize STC file time start
    stc(ii).tstep = 1/sfreq;                                                %initialize STC file sampling
end
stc(1).vertices = [];   stc(2).vertices = [];                               %Initialize for Simulation
t = tot_t_time/1000*sfreq;                                                  %# time points
time = (1/sfreq:1/sfreq:tot_t_time/1000); disp(['length of time vector...', num2str(length(time))]); %time vector [sec]

% Initialize Fwd Solutions
scalarX = double(fwd.nsource/fwd4.nsource);                                 %ratio for hypersampled and ico-4 amplitudes
G = [];                                                                     %Initialize for simulation
Gvpl_orig = sdiv_fwd(vpl_sdiv).rawG;                                        %equals sdiv_fwd(80).rawG in inv folder
if ~exist('raw_vpl_svd.mat','file')
    [~,Ur,Sr,Vr,~,~,~] = nmra(Gvpl_orig,ploton,savemem,p_sing,nthresh);     %svds of raw forward solution
    save('raw_vpl_svd.mat','Ur','Sr','Vr');                                 %save SVD outputs so no redo required
else
    load('raw_vpl_svd.mat');                                                %load SVD outputs    
end
Gvpl = Ur(:,1)*Sr(1,1);                                                     %reduce fwd solution to most significant mode

% Initialize Channel Nos and Noise Covariance
coil_type = [meas.info.chs(:).coil_type];
magnetom = ismember(coil_type,3022);    magn = find(magnetom==1);   N_m = sum(magnetom); 
gradiom = ismember(coil_type,3012);     grad = find(gradiom==1);    N_g = sum(gradiom);  
C0 = zeros(N_m + N_g, N_m + N_g);                                   

%% SIMULATE CORTICAL COMPONENTS OF MEDIAN-NERVE ELECTRICAL STIMULATION ERP
% Specify Cortical Currents with Gabor Atoms
gabor_freq = 1000;
gabor(1,:) = Gatom(0.30,0,40,gabor_freq,12,cort_t_time, cort_repeat_delay) + ...%cS1
Gatom(-1.06,0,40,gabor_freq,25,cort_t_time, cort_repeat_delay) + Gatom(-0.80,0,400,gabor_freq,75,cort_t_time, cort_repeat_delay)...
+ Gatom(-0.25,0,800,gabor_freq,160,cort_t_time, cort_repeat_delay) + Gatom(0.12,0,400,gabor_freq,212,cort_t_time, cort_repeat_delay); 
gabor(2,:) = Gatom(0.03,0,500,gabor_freq,0,cort_t_time, cort_repeat_delay) + ...%cS2
Gatom(0.40,0,400,gabor_freq,75,cort_t_time, cort_repeat_delay) + Gatom(-0.30,0,600,gabor_freq,150,cort_t_time, cort_repeat_delay) + ...
Gatom(0.12,0,800,gabor_freq,250,cort_t_time, cort_repeat_delay);
gabor(3,:) = Gatom(0.02,0,400,gabor_freq,0,cort_t_time, cort_repeat_delay) + ...%PPC
Gatom(0.24,0,80,gabor_freq,75,cort_t_time, cort_repeat_delay) + Gatom(0.6,0,200,gabor_freq,160,cort_t_time, cort_repeat_delay) + ...
Gatom(0.1,0,200,gabor_freq,225,cort_t_time, cort_repeat_delay);
gabor(4,:) = Gatom(0.12,0,200,gabor_freq,75,cort_t_time, cort_repeat_delay) + ...%iS2
Gatom(0.25,0,800,gabor_freq,125,cort_t_time, cort_repeat_delay) + Gatom(-0.25,0,400,gabor_freq,175,cort_t_time, cort_repeat_delay);

% Simulate Gabor-Based Activity in Each Region of Interest
nuse = length(fwd.src(1).vertno);
for ii = 1:length(label_string)
    if ii <= 3                                                              %cS1, cS2 and PPC  
        label = mne_read_label_file(fullfile(strcat(datapath,'/label/',prefix,'-ico-2/lh/',label_string{ii})));
        [vertno,~,X(ii).lk] = intersect(label.vertices+1,fwd.src(1).vertno);%keep track of dipole indices
        stc(1).vertices = cat(2,stc(1).vertices,vertno);                    %assign STC vertices    
        dist = 1000*src(1).nearest_dist(vertno);                            %for distance decaying profile to dipoles
    else                                                                    %iS2
        label = mne_read_label_file(fullfile(strcat(datapath,'/label/',prefix,'-ico-2/rh/',label_string{ii})));
        [vertno2,~,ind] = intersect(label.vertices+1,fwd.src(2).vertno); X(ii).lk = ind + nuse; %track dipole indices
        stc(2).vertices = cat(2,stc(2).vertices,vertno2);                   %assign STC vertices
        dist = 1000*src(2).nearest_dist(vertno2);                           %for distance decaying profile to dipoles
    end
    
    % Pick out forward vectors corresponding to dipoles
    G = cat(2,G,fwd.sol.data(:,X(ii).lk));                                  %raw forward solutions
    X(ii).data = zeros(length(X(ii).lk), t);                                %initialize X(ii) data matrix
    
    % Apply distance decaying profile to dipole activity
    cst_ind = find(time == scort_t_time/1000) + 1;                          %time to start cortical currents
    cen_ind = t;
    wts = repmat((1-(dist/max(dist)).^2/2)',1,length(gabor(ii,:)));         %distance decaying weights
    gabor_space = 10*10^-9*wts.*repmat(gabor(ii,:),size(wts,1),1);          %resample gabor profiles in space
    X(ii).data(:,cst_ind:cen_ind) = resample(gabor_space',sfreq/gabor_freq,1)'; %resample gabor profiles in time
    X(ii).data = X(ii).data(:,1:end-trim_end);                              %trim out last part as it has weird behavior
    
    % Pick out orientation vectors corresponding to dipoles
    X(ii).orient = fwd.source_nn(X(ii).lk,:);
end
Xhyper = cat(1,X(:).data);                                                  %simulated activity in hypersampled dipoles

%% SIMULATE SUBCORTICAL COMPONENT OF MEDIAN-NERVE ELECTRICAL STIMULATION ERP
% Specify Subcortical Current: Many Overlap
num_scurr = 1.5; f = 100; phi = pi/3;                                       %freq, phase shift for cos^2 function
Xvpl = zeros(1,size(Xhyper,2));                                             %initialize subcortical currents
scurr_str = th_cd*sdiv_fwd(vpl_sdiv).volume; amp = num_scurr*scurr_str;     %put strength of all dipoles into this mode
t_ind = scort_t_time/1000*sfreq;                                            %time range for 1st cos^2 function
mshape = amp*cos(2*pi*f*time(1:t_ind) + phi).^2;                            %alternative: sinc(linspace(-5,10,delay))
Xvpl(7:37) = mshape(3:33);                                                  %single series
sim_on{1} = [7 37];							    %simulation on points
if num_scort_repeats > 1
    t_rep = scort_repeat_delay/1000*sfreq;                                  %repeat range for cos^2 function
    for i = 1:num_scort_repeats-1
        Xvpl(7+t_rep*i-1:37+t_rep*i-1) = mshape(3:33);                      %repeated series
        sim_on{i+1} = [7+t_rep*i-1, 37+t_rep*i-1];                          %simulation on points
    end   
end
Xvpl = smooth(Xvpl)';                                                       %smooth out

%% COMBINE CORTICAL, SUBCORTICAL, AND NOISE SIMULATIONS
Ynoisefree = Gvpl*Xvpl + G*Xhyper;                                          %Meas from subcortical + cortical currents      
cort_t_time = cort_t_time - trim_end/sfreq*1000;                            %update duration of cortical gabors [msec]
tot_t_time = scort_t_time + cort_t_time;                                    %update duration of full simulation [msec]
t = tot_t_time/1000*sfreq;                                                  %update # time points
time = time(1:end-trim_end);                                                %update time vector [sec]

%% SIMULATE NOISE SERIES AND COVARIANCE MATRIX
% White Gaussian noise based on fixed SNR value and signal amplitude
sigmaNu_m = norm(Ynoisefree(magnetom,:),'fro')^2/(N_m*t*SNR);               %Magnetometer noise variance
sigmaNu_g = norm(Ynoisefree(gradiom,:),'fro')^2/(N_g*t*SNR);                %Gradiometer noise variance
nu_m = normrnd(0, sqrt(sigmaNu_m), N_m, t);                                 %Magnetometer noise series
nu_g = normrnd(0, sqrt(sigmaNu_g), N_g, t);                                 %Gradiometer noise series

% Add Noise to MEG Observations
Y(magnetom,:) = Ynoisefree(magnetom,:) + nu_m;                              %Noisy magnetometer Measurement
Y(gradiom,:) = Ynoisefree(gradiom,:) + nu_g;                                %Noisy gradiometer Measurement

% Noise Covariance Matrix
for i = 1:N_m
    C0(magn(i), magn(i)) = sigmaNu_m;                                       %Magnetometer noise covariance
end
for i = 1:N_g
    C0(grad(i), grad(i)) = sigmaNu_g;                                       %Gradiometer noise covariance
end

% Scale Measurements and Currents
scalarY = max(max(abs(meas.evoked.epochs(1:306,t0:t0+300))))/max(max(abs(Y)));%scale factor  for realistic amplitude level
Ynoisefree = scalarY*Ynoisefree;                                            %scale down noisefree measurement
Y = scalarY*Y;                                                              %Scale down noisy measurement
C0 = scalarY^2*C0;                                                          %For Consistency: as noise level goes down automatically
for ii = 1:length(label_string)
    X(ii).data = scalarY*X(ii).data;                                        %For Consistency: Est. Xs will be scaled by scalarY                      
end
Xhyper = scalarY*Xhyper;                                                    %For Consistency: Est. Xs will be scaled by scalarY
Xvpl = scalarY*Xvpl;                                                        %For Consistency: Est. Xs will be scaled by scalarY

%% COMPUTE RESULTANTS AND STORE FOR PLOTTING
% Put All Subcortical Currents in Hypersampled Dipole Space
Xvpl_plot.numdipoles = sdiv_fwd(vpl_sdiv).ndipoles;                         %VPL # dipoles - also size(Vr,1)/3; 
Xvpl_plot.mode = Xvpl;                                                      %VPL current in mode space (nAm)
sdip = Vr(:,1)*Xvpl;                                                        %project subcortical mode to (hypersampled) dipole space 
Xvpl_plot.dip_or(:,:,1) = sdip(1:3:end,:);                                  %orientation 1
Xvpl_plot.dip_or(:,:,2) = sdip(2:3:end,:);                                  %orientation 2
Xvpl_plot.dip_or(:,:,3) = sdip(3:3:end,:);                                  %orientation 3
Xvpl_plot.dip = sqrt(sum(Xvpl_plot.dip_or.^2,3));                           %dipole current estimate (amplitude)
[Xvpl_plot.dip_res, Xvpl_plot.dip_resor] = vec_sum(Xvpl_plot.dip_or);       %net activity per dipole

% Put All Cortical Currents in Hypersampled Dipole Space
Xc_plot(1:length(label_string)) = struct('numdipoles',[],'cstrdip',[],'dip',[],'orient',[],'dip_or',[],'dip_res',[],'dip_resor',[]);
Xc_hyper_plot(1:length(label_string)) = struct('dip',[],'res',[]); 
for ii = 1:length(label_string)
    Xc_plot(ii).numdipoles = size(X(ii).data,1);                            %Cortical Region # Dipoles
    Xc_plot(ii).dip = X(ii).data;                                           %cortical amplitudes in hypersampled space
    Xc_plot(ii).orient = X(ii).orient;                                      %orientations of dipoles
    nn(:,:,1) = Xc_plot(ii).orient; nn = repmat(nn,[1,1,t]);                %orientations    
    Xc_plot(ii).dip_or(:,:,1) = squeeze(nn(:,1,:)).*Xc_plot(ii).dip;        %1 X #timepoints
    Xc_plot(ii).dip_or(:,:,2) = squeeze(nn(:,2,:)).*Xc_plot(ii).dip;        %1 X #timepoints
    Xc_plot(ii).dip_or(:,:,3) = squeeze(nn(:,3,:)).*Xc_plot(ii).dip;        %1 X #timepoints
    [Xc_plot(ii).dip_res, Xc_plot(ii).dip_resor] = vec_sum(Xc_plot(ii).dip_or);     %net activity (amplitude and orientation) across all dipoles
    Xc_hyper_plot(ii).res = repmat(Xc_plot(ii).dip_res,Xc_plot(ii).numdipoles,1);   %net activity across all dipoles repeated
    clear nn
end
Xhyper_plot.dip = cat(1,Xc_plot(:).dip);                                    %dipole currents           
Xhyper_plot.dip_res = cat(1,Xc_hyper_plot(:).res);                          %resultant currents

%% FINAL SAVE AND PLOT
% Write STC Files (hypersampled source space) of Simulated Cortical Activity
X_lh = cat(1,X(1:3).data); M = size(X_lh,1);
stc(1).data = Xhyper_plot.dip_res(1:M,:);
stc(2).data = Xhyper_plot.dip_res(M+1:end,:);
mne_write_stc_file1(strcat(savename,'-lh.stc'),stc(1))
mne_write_stc_file1(strcat(savename,'-rh.stc'),stc(2))

% Save Time Series and Simulation Parameters
if saveon
    save(strcat(savename,'.mat'),...
    'N','time','meas','magnetom','gradiom',...
    'sfreq', 'scort_t_time','cort_t_time','tot_t_time', 'cort_repeat_delay','t', ...
    'SNR','p_sing', 'nthresh','t0','vpl_sdiv', 'scalarX','scalarY',...
    'G','Ur','Sr','Vr','Gvpl',...
    'X','Xhyper','Xvpl','sim_on',...
    'C0','Y','Ynoisefree',...
    'Xvpl_plot','Xc_plot','Xhyper_plot','stc', '-v7.3')
end
disp('completed saving simulation');

%Plot Time Series in MATLAB
%figure, set(gcf,'color','white'); hold all;
%plot(time*1000, (10^9)*Xhyper_plot.dip_res'); plot(time*1000, (10^9)*Xvpl_plot.dip_res','k');
%ylabel('Source Current [nAm]'); xlabel('Time (Milliseconds)'); 
diary off;

ploton = 0;                                                                 % plot covs (1) and redo svd (1)
SM04_simSEP_Cw(prefix, meastype, datatype, ploton, overlap);                % compute whitener etc and append to measfile
SM04_simplotfinal(overlap);                                                 % plot simulation results
sensor_space_visualization;						    % plot sensor space field map