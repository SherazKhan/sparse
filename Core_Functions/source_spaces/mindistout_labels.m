function mindistout_labels(prefix, fwd, fwd_mindist, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
% Edited by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ID Vertices that are in Fwd but not in Mindist Fwd
%l1 = min(length(fwd.src(1).vertno), length(fwd_mindist.src(1).inuse));
lk = fwd.src(1).vertno(~fwd_mindist.src(1).inuse);%(1:l1));
rows = length(lk);
%l2 = min(length(fwd.src(2).vertno), length(fwd_mindist.src(2).inuse));
lk2 = fwd.src(2).vertno(~fwd_mindist.src(2).inuse);%(1:l2));
rows2 = length(lk2);

%% Store Labels that Fail Left Mindist Requirement
if nargin == 3
    fid = fopen(strcat(prefix,'-meg-all-mindistout-lh-label.txt'),'w');
else
    fid = fopen(varargin{1},'w');
end
for ii = 1:rows
    if ii < rows
        fprintf(fid,[num2str(lk(ii)) '\n']);
    else
        fprintf(fid,num2str(lk(ii)));
    end
end
fclose(fid);

%% Store Labels that Fail Right Mindist Requirement
if nargin == 3
    fid = fopen(strcat(prefix,'-meg-all-mindistout-rh-label.txt'),'w');
else
    fid = fopen(varargin{2},'w');
end
for ii = 1:rows2
    if ii < rows2
        fprintf(fid,[num2str(lk2(ii)) '\n']);
    else
        fprintf(fid,num2str(lk2(ii)));
    end
end
fclose(fid);

end
