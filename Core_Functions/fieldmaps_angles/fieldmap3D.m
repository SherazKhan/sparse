function [y,sel] = fieldmap3D(data, fwd, fieldmap_fname, sel, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a field map with a given fwd solution and source current distn
% Done by modifying file name specified by data into fieldmap_fname
% Written by Pavitra Krishnaswamy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Source Currents
if nargin == 4
    x = ones(size(fwd,2),1); 
else
    x = varargin{1};
end

% Read in Data File to Modify
newdata.info = data.info;                               % copy information and channel fields
if length(data.evoked) > 1
    newdata.evoked = data.evoked(1);                    % copy evoked fields for first response only
else
    newdata.evoked = data.evoked;                       % copy evoked fields
end

% Generate No-Noise Measurement from Cortical Lead Field
t = newdata.evoked.times;                               % ERP recording times
y = fwd*x;                                              % measurement
if size(fwd,1) ~= length(sel)
    y = y(sel,:);                                       % only include y for selected channels
end
newdata.evoked.epochs = repmat(y, 1, length(t));        % repeat across time

% Write Modified Evoked File
fiff_write_evoked(fieldmap_fname, newdata);             % write field map file

end