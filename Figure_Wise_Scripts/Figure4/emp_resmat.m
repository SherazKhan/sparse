function K = emp_resmat(xc, xs, lk, clk, cp, regname, ssubdiv, svolume, Gtheta, trimP_best, X_best, ploton)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Empirical Resolution Matrix
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simulation Parameters and Current Quantum 
clk_inds = 1:1:length(clk);                                                 % Indices of cortical elements
num_cort = length(clk);                                                     % Number of cortical elements
slk_inds  = num_cort + 1:length(trimP_best);                                % Indices of subcortical elements
xsim(clk_inds) = xc; xsim(slk_inds) = xs(slk_inds);                         % Simulated Currents, Already Scaled by CurrStr [Am]

%% Extract Location of Estimated Current
for j = 1:length(trimP_best)
  if trimP_best{j} <= num_cort                                              % Estimated location is cortical    
    est(j) = ceil((find(trimP_best{j} == clk_inds))/cp);                    % Estimated Patch Index
  else                                                                      % Estimated location is subcortical
    count = 0;                                                              % Start Counter                                        
    for sreg = 1:length(regname)    
      for sdiv = 1:length(ssubdiv{sreg})
        count = count + 1;                                                  % Counter for Subcortical Divisions
    	inds = find(svolume{sreg,sdiv} == lk(trimP_best{j}));               % ID anatomical region, subdivision index
        if isempty(inds) == 0                                               % inds not empty, i.e. this is the region
        est(j) = count + num_cort/cp;                                       % count is division #, add to cort patch index to derive matrix element                             
        end
      end
    end
  end
end

%% Compute Resolution Matrix
% Pavitra Initial Approach
% Kp = zeros(length(est)); 
% for i = 1:length(est)
%     a = zeros(1,length(est));                                             % initialize ith column of resolution matrix
%     x_curr = X_best{i}{end};                                              % Estimated Current in Division i, Scaled by CurrStr
%     a(est(i)) = x_curr(x_curr ~= 0)/xsim{i};                              % the value of this column is 1 if equal, else not
%     Kp(:,i) = a;                                                          % Only 1 non-zero value for resolution matrix
% end

% Mattis Suggested Formula: Does not Explicitly Use x_curr Output by SP, Rather Computes it with Least Squares
K = zeros(length(est)); 
for j = 1: length(est)                                                      % Index of Actual Location
    k = est(j);                                                             % Index of Estimated Location                                                  
    Knum = Gtheta(:,k)'*Gtheta(:,j);                                        % Numerator for Least Squares Fit
    Kdenom = Gtheta(:,k)'*Gtheta(:,k);                                      % Denominator for Least Squares Fit
    K(k,j) = Knum/Kdenom;
end

%% Plot Resolution Matrix
if ploton
  figure, set(gcf,'color','white'); 
  imagesc(abs(K)); colormap(flipud('gray'));                                % Plot with MATLAB standard colorbar and no embellishment
  %imagesc(abs(Kp)); colormap(flipud('gray'));                              % Plot with MATLAB standard colorbar and no embellishment
  %isequal(K, Kp)
end

end