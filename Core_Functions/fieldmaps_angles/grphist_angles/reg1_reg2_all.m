function angf1f2 = reg1_reg2_all(reg1, reg2, flag)

%% Collect Eigenmodes for all Divisions in Region 1
U1full = []; S1full = []; 
for j = 1:length(reg1)
    if ~isempty(reg1(j).U)
      if strcmp(flag,'all')
        U1full = [U1full reg1(j).U];
        %S1full = [S1full reg1(j).S]; %given different #eigenmodes difficult
      elseif strcmp(flag,'rand')
        mode_sel1 = randperm(size(reg1(j).U,2), 1);
        U1full = [U1full reg1(j).U(:,mode_sel1)];
      elseif strcmp(flag,'top')
        U1full = [U1full reg1(j).U(:,1)];
      end
    end
end

%% Collect Eigenmodes for all Divisions in Region 2
U2full = []; S2full = []; 
for j = 1:length(reg2)
    if ~isempty(reg2(j).U)
	  if strcmp(flag,'all')
        U2full = [U2full reg2(j).U];
        %S2full = [S2full reg2(j).S];
      elseif strcmp(flag,'rand')
        mode_sel2 = randperm(size(reg2(j).U,2), 1);
        U2full = [U2full reg2(j).U(:,mode_sel2)];
      elseif strcmp(flag,'top')
        U2full = [U2full reg2(j).U(:,1)];
      end
    end
end

%% Angles Between Gain Matrices of Regions 1 and 2
if isequal(reg1, reg2)
  angf1f2 = [];
else
  %rows to exclude from U2full and U1full to make mutually exclusive
  [~, excl_rows2, excl_rows1] = intersect(U2full', U1full','rows');   	
  incl_rows2 = setdiff([1:1:size(U2full,1)], excl_rows2); 
  incl_rows1 = setdiff([1:1:size(U1full,1)], excl_rows1); 
  angf1f2 = subspacea(U2full(incl_rows2,:), U1full(incl_rows1,:))*180/pi;
end
clear U1full U2full %S1full S2full

end