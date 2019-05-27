function [Gtheta, Ur, Sr, Vr, best_numeigen, nmr, d_nmra] = nmra(G,ploton,savemem,p_sing,th)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized Mean Representation Accuracy (Limpiti 2006) for Forward Matrices
% Savemem is a flag for when we want to save memory 
% If savemem - compute reduced svd with 35 largest modes only and nmra = trace(S'*S)/trace(S'*S)
% No savemem - compute full svd and set nmra = trace(G'*Ufull*Ufull'*G)/trac(G'*G)
% Determine best # modes by summing nmra's till reach threshold MRA
% Take reduced svd with 10 modes to return Gtheta, Ur, Sr, Vr for analysis
% Written by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu) 
% With Inputs from Gabriel Obregon (obregon@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FORM DISJOINT LEAD FIELDS
if savemem == 0
    [Ufull,Sfull,Vfull] = svd(full(G)); 
    N = size(Ufull,2);
else
    S = svds(G,35,'L');     %for memory convenience
    N = min(35, length(S)); 
end

%% COMPUTE NORMALIZED MEAN REPRESENTATION ACCURACY
disp(['Computing NMRA using p = 1:',num2str(p_sing),' modes...'])
for p = 1:N
    if savemem == 0
        nmr(p) = trace(G'*Ufull(:,1:p)*Ufull(:,1:p)'*G)/trace(G'*G);
    else
        nmr(p) = trace(S(1:p)'*S(1:p))/trace(S'*S);
    end
    cumsum(p) = sum(nmr(1:p));
end

%% PLOT  NMRA
if ploton
    figure, plot(nmr','x'),grid on;
    figure, plot(cumsum,'o'); grid on;
    %disp(cumsum([1:5]))
end

%% CHOOSE BEST NUMBER OF EIGENMODES
%d = diff([0 cumsum]); %disp(d([1:5]))
d_nmra = diff([0 nmr])*100; %disp(d_nmra(1:10))
best_numeigen = find(nmr > th/100,1);

%% RETURN REDUCED SVD FOR TOP p_sing MODES 
if savemem == 0
    numeigen = min(p_sing,size(Sfull,2)); 
    Ur = Ufull(:,1:numeigen);
    Sr = Sfull(1:numeigen,1:numeigen);
    Vr = Vfull(:,1:numeigen);
else
    numeigen = min(p_sing,length(S)); 
    [Ur,Sr,Vr] = svds(G,numeigen,'L');
end
Gtheta = Ur*Sr;
    
end 