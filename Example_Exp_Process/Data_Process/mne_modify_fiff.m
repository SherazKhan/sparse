function mne_modify_fiff(infile, outfile, newdata)
% Adapted from mne_ex_read_write_raw by Matti Hamalainen, MGH Martinos Center
% License : BSD 3-clause
% Read in Data from Fiff File and Modify to NewData
% Modified by Pavitra Krishnasawmy - 04/2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end
me = 'MNE:mne_modify_fiff';
if nargin ~= 3
    error(me,'Incorrect number of arguments');
end

%% Setup for Reading the Raw Data
try
    raw = fiff_setup_read_raw(infile);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end

%% Setup Pick List
picks = 1:1:size(newdata,1);
[outfid,cals] = fiff_start_writing_raw(outfile,raw.info,picks);

%% Setup the Reading Parameters
from        = raw.first_samp;
to          = raw.last_samp;
quantum_sec = 0.3; quantum     = ceil(quantum_sec*raw.info.sfreq);
%quantum     = to - from + 1;        %to read the whole file at once

%% Read and write all the data
first_buffer = true;
for first = from:quantum:to
    last = first+quantum-1;
    if last > to
        last = to;
    end
    try
        [data, times] = fiff_read_raw_segment(raw,first,last,picks);
    catch
        fclose(raw.fid);
        fclose(outfid);
        error(me,'%s',mne_omit_first_line(lasterr));
    end
    
    %% Modify Data to NewData
    data = newdata(picks,first-from+1:last-from+1);
    %if newdata same as Y or data - checked the below
    %isequal(data, newdata(picks,first-from+1:last-from+1)) is 1
    
    %% Writing
    fprintf(1,'Writing...');
    if first_buffer
       if first > 0
           fiff_write_int(outfid,FIFF.FIFF_FIRST_SAMPLE,first);
       end
       first_buffer = false;
    end
    fiff_write_raw_buffer(outfid,data,cals);
    fprintf(1,'[done]\n');
end

fiff_finish_writing_raw(outfid);
fclose('all'); %fclose(raw.fid);
end