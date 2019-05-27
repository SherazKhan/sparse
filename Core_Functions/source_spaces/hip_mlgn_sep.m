clear; clc; close all;
dpath = '/autofs/eris/purdongp/users/pavitrak/sourceloc_data/AEP/';
subject = 'manny';
addpath(strcat(dpath,subject,'/src'));
saveon = 1; 

lfull = importdata(strcat(subject,'_lh.txt')); 
lsz = size(lfull)
lmgn = importdata(strcat(subject,'_lmgn.txt'));
llgn = importdata(strcat(subject,'_llgn.txt'));

rfull = importdata(strcat(subject,'_rh.txt'));
rsz = size(rfull)
rmgn = importdata(strcat(subject,'_rmgn.txt'));
rlgn = importdata(strcat(subject,'_rlgn.txt'));

i_lmgn = ismember(lfull,lmgn,'rows'); 
i_llgn = ismember(lfull,llgn,'rows'); 

i_rmgn = ismember(rfull,rmgn,'rows'); 
i_rlgn = ismember(rfull,rlgn,'rows'); 

lfull = lfull(~i_lmgn & ~i_llgn,:);
lsz = size(lfull)
rfull = rfull(~i_rmgn & ~i_rlgn,:);
rsz = size(rfull)

savename = strcat(dpath,subject,'/src/',subject,'_lh.txt');
if saveon
    save(savename, 'lfull', '-ASCII'); %type(savename);
end

savename = strcat(dpath,subject,'/src/',subject,'_rh.txt');
if saveon
    save(savename, 'rfull', '-ASCII'); %type(savename);
end

load(strcat(dpath, subject, '/src/',subject,'_lh_geom.mat'),'volume','ndipoles');
disp('left'); disp(volume); disp(ndipoles);
volume = (lsz(1)/3)/ndipoles*volume; ndipoles = lsz(1)/3;
if saveon
    save(strcat(dpath, subject, '/src/',subject,'_lh_geom.mat'),'volume','ndipoles');
end
clear volume ndipoles

load(strcat(dpath, subject, '/src/',subject,'_rh_geom.mat'),'volume','ndipoles');
disp('right'); disp(volume); disp(ndipoles);
volume = (rsz(1)/3)/ndipoles*volume; ndipoles = rsz(1)/3;
if saveon
    save(strcat(dpath, subject, '/src/',subject,'_rh_geom.mat'),'volume','ndipoles');
end
clear volume ndipoles lfull rfull lsz rsz i_lmgn i_rmgn i_rlgn i_llgn