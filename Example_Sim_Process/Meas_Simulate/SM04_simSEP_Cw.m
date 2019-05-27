function SM04_simSEP_Cw(prefix, meastype, datatype, ploton,  overlap)

%% File Names
datapath = strcat('/autofs/eris/purdongp/users/pavitrak/sourceloc_data/',datatype,'/',prefix,'/');
meas_fname = strcat(datapath, 'meas/', prefix,'_vpl_somato_meg_sim',overlap,'.mat');  %meas file
fwd_fname = strcat(datapath, 'fwd/', prefix,'-meg-all-fixed-fwd.fif');      %fwd file - uses 3-BEM, no mindist
diary(strcat(prefix,'_',meastype,'_sim_whitener_matlog',overlap));          %log of simulation data

%% Make Channel Selectors
load(meas_fname, 'C0','Y','stc');                                           % load cov, Y and stc structures
sel = 1:1:306;                                                              % all channels good for simulated data 
Y = Y(sel,:);                                                               % select only good channels

%% Measurement (Simulated)
N = size(Y,1);                                                              % # channels
t = size(Y,2); t1 = 1; t2 = t;                                              % # timepoints
time = linspace(0,t-1,t)*stc(1).tstep;                                      % time vector
%time = linspace(ceil(Fs*stc(1).tmin),ceil(Fs*stc(1).tmin)+(t-1),t)*stc(1).tstep;

%% Noise Covariance (Simulated)
Leff = 1;                                                                   % number of epochs
P = eye(N);                                                                 % projector is identity
Cn = C0/Leff;                                                               % simulated covariance
Cp = P*(Cn)*P;                                                              % projected covariance
fwd = mne_read_forward_solution(fwd_fname);                                 %cortical forward file
coil_type = [fwd.chs(:).coil_type];                                         % obtain coil type information
cij = ismember(coil_type,3024);                                             % demarcate magnetometers

%% Load or Regularize Noise Covariance
Cl = Cp;                                                                    % loaded covariance

%% Derive EigenDecomposition of Noise Covariance and Compute Whitener 
diag_std = sqrt(diag(Cl));                                                  % noise stdevs ordered correctly
Cw = diag(1./diag_std);                                                     % whitener - preserves ordering as diagonal
Cw = double(Cw);                                                            % need double precision - svds for later

%% Plot Covs and Eigenvalues
if ploton == 1
    figure, set(gcf,'color','white'); imagesc(Cl(cij,cij)); colorbar;
    figure, set(gcf,'color','white'); imagesc(Cl(~cij,~cij)); colorbar;
    figure, set(gcf,'color','white'); semilogy(sort(real(lambda2p(lambda2p>lamthresh)))); 
    %Rc = chol(Cl); Cw = inv(Rc);  %previous case
    pause; close;
end

%% Whiten Data 
Yw = Cw*Y;                                                                  % whitened measurement 
C = eye(N);                                                                 % noise cov after whitening

%% Save Important Variables and Call Inverse Solutions
savename = meas_fname; 
save(savename, 'prefix','meastype','datatype',...
               'N','Leff','cij','sel','t1','t2','t','stc','time','Yw',...
               'C0','Cn','P','Cp','Cl','Cw','C','-append');    
diary off

end