function [Xls_scaled,Xls,logLs,Rdiag] = source_estimates(alg,C,G,N,Y,t1,t2,SNR,stop,currstr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
% Edited by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute Source Covariance -> SNR of Whitened Data is ~ 1/sigmaX 
sigmaX = N/trace(G*G');
R0 = sigmaX*eye(size(G,2));
t = t2-t1+1;

%% Set Regularization of Parameter
lambda_sq = 1/SNR;

%% Run Inverse Solution with Chosen Algorithm
switch alg 
    case 'MNE'
        % Form whitened and regularized inverse operator
        W = (R0*G')/(G*R0*G'+lambda_sq*C);
        % Compute MNE source estimates
        Xls = W*Y;
    case 'dSPM'
        % Form inverse operator
        W = (R0*G')/(G*R0*G'+lambda_sq*C);
        % Compute dSPM source estimates
        Xls = (W*Y)./repmat(sqrt(diag(W*C*W')),1,t);
    case 'sMAP-EM'
        % Compute sMAP-EM source estimates
        [Xls,logLs,Rdiag] = sMAP_EM(C,G,lambda_sq, R0,stop,t1,t2,Y);
    otherwise
        error('Unknown Source Localization Method')
end

%% Scale Current Estimates by Current Strengths
for i = 1:size(Xls,1)
    Xls_scaled(i,:) = Xls(i,:)*currstr(i);
end

end