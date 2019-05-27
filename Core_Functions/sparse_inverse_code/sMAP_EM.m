function [Xls,logLs,Rdiag] = sMAP_EM(C,G,lambda_sq,R0,stop,t1,t2,Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATIC MAXIMUM A POSTERIORI EXPECTATION MAXIMIZATION (sMAP-EM) ALGORITHM
% Written by Camilo Lamus (lamus@mit.edu)
% Recent edits by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
% Edit by Pavitra Krishnaswamy (pavitrak@mit.edu) - 1:T replaced by t1:t2  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
R = R0/lambda_sq;
Rdiag = diag(R);
T = t2-t1+1;

%% Set Prior
% Prior on the variance of tap input noise of dipole sources (R(m,m)) is set to inverse gamma with parameters alpha and beta
% The mean on the prior is beta/(alpha-1), the value of the variance of the sources in MNE
o = R(1,1);
l = 1/o;
prior.alpha = o^2/l+2;
prior.beta = o*(o^2/l+1);

%% STATIC EM-MAP LOOP
% obtain MNE solutions 
% find maximum in conditional expected value(log of the joint complete data and parameters)
% given the observed data and parameters in previous iteration
% Inversion using MNE and computation of MSE
for ii = 1:stop
    % iterationEM_static = ii
    
    % Expectation (MNE)
    covY = G*R*G' + C;
    W = (R*G')/(covY);
    Xls = W*Y;
        
    % Minimization
    P = (R-W*G*R);
    B = T*P;
    logLs(ii) = -T/2*log(det(covY)); %logLs(ii) = -T/2*find_log_det(covY);
    for t = t1:t2
        B = B + Xls(:,t)*Xls(:,t)';
        logLs(ii) = logLs(ii) - (Y(:,t)'/covY)*Y(:,t)/2;
    end
    Rdiag(:,ii+1) = (diag(B)/2+prior.beta)/(T/2+(prior.alpha+1));
    R = diag(Rdiag(:,ii+1));
end

end