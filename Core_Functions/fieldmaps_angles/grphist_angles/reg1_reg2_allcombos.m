function [angf1f2_allcombos, varargout] = reg1_reg2_allcombos(reg1, reg2, type_flag, reg1p, reg2p)

%% Collect Eigenmodes for all Region 1 Subdivisions
U1full = []; S1full = []; 
for j = 1:length(reg1)
    if ~isempty(reg1(j).U)
      if strcmp(type_flag,'all') || strcmp(type_flag,'all_top')
        U1full = [U1full reg1(j).U];
        S1full = [S1full diag(reg1(j).S)'];
      elseif strcmp(type_flag,'rand')
        mode_sel1 = randperm(size(reg1(j).U,2), 1);
        U1full = [U1full reg1(j).U(:,mode_sel1)];
        S1full = [S1full reg1(j).S(mode_sel1, mode_sel1)'];
      elseif strcmp(type_flag,'top')
        U1full = [U1full reg1(j).U(:,1)];
        S1full = [S1full reg1(j).S(1,1)'];
      end
    end
end

%% List all Possible Subsets of Region 1 Modes
%ignore empty subsets
U1_subs = []; U1_subs_inds = [];
sz1 = size(U1full,2);
mode_inds1 = 1:1:sz1;
for k1 = 1:length(mode_inds1)
    U1_subs{k1} = nchoosek(mode_inds1,k1);
end
%reorganize
count = 1;
for kk = 1:length(U1_subs)
    for ii = 1:size(U1_subs{kk},1)
        U1_subs_inds{count} = U1_subs{kk}(ii,:);
        count = count+1;
    end
end
if strcmp(type_flag, 'all_top')
    for jj = 1:reg1p
        t_mode_inds1(jj,:) = jj:reg1p:sz1;
    end
    [U1_subs_inds, flag_red1] = reduce_combos(U1_subs_inds, t_mode_inds1);
end

%% Collect Eigenmodes for all Region 2 Subdivisions
U2full = []; S2full = []; 
for j = 1:length(reg2)
    if ~isempty(reg2(j).U)
      if strcmp(type_flag,'all') || strcmp(type_flag,'all_top')
        U2full = [U2full reg2(j).U];
        S2full = [S2full diag(reg2(j).S)'];
      elseif strcmp(type_flag,'rand')
        mode_sel2 = randperm(size(reg2(j).U,2), 1);
        U2full = [U2full reg2(j).U(:,mode_sel2)];
        S2full = [S2full reg2(j).S(mode_sel2,mode_sel2)'];
      elseif strcmp(type_flag,'top')
        U2full = [U2full reg2(j).U(:,1)];
        S2full = [S2full reg2(j).S(1,1)'];
      end
    end
end

%% List all Possible Subsets of Region 2 Modes
%ignore empty subsets
U2_subs = []; U2_subs_inds = [];
sz2 = size(U2full,2);
mode_inds2 = 1:1:sz2;
for k2 = 1:length(mode_inds2)
    U2_subs{k2} = nchoosek(mode_inds2,k2);
end
%reorganize
count = 1;
for kk = 1:length(U2_subs)
    for ii = 1:size(U2_subs{kk},1)
        U2_subs_inds{count} = U2_subs{kk}(ii,:);
        count = count+1;
    end
end
if strcmp(type_flag, 'all_top')
    for jj = 1:reg2p
        t_mode_inds2(jj,:) = jj:reg2p:sz2;
    end
    [U2_subs_inds, flag_red2] = reduce_combos(U2_subs_inds, t_mode_inds2);
end

%% Angles Between Each Combination of Region 1 and Region 2 Modes
angf1f2 = []; angf1f2_allcombos = [];
for k2 = 1:length(U2_subs_inds)
    for k1 = 1:length(U1_subs_inds)
        cols_select2 = U2_subs_inds{k2};
        cols_select1 = U1_subs_inds{k1};
        U2_subset = U2full(:,cols_select2);
        U1_subset = U1full(:,cols_select1);
        if isequal(reg1, reg2)
          if isempty(intersect(cols_select2, cols_select1))                 %mutually exclusive
          %if ~isequal(U2_subset, U1_subset)                                %not identical
            angf1f2{k2, k1} = subspacea(U2_subset, U1_subset)*180/pi;
            angf1f2_allcombos = [angf1f2_allcombos angf1f2{k2, k1}(:)'];
          end
        else
            angf1f2{k2, k1} = subspacea(U2_subset, U1_subset)*180/pi;
            angf1f2_allcombos = [angf1f2_allcombos angf1f2{k2, k1}(:)'];
        end
    end
end

%% Output Details and Clear
varargout{1} = angf1f2; 
varargout{2} = U1_subs_inds;    varargout{3} = U1full;  varargout{4} = S1full;
varargout{5} = U2_subs_inds;    varargout{6} = U2full;  varargout{7} = S2full;
%save(strcat('test_',type_flag,'.mat'),'angf1f2','U1_subs_inds','U2_subs_inds','U1full','U2full')
clear U1full U2full S1full S2full

end