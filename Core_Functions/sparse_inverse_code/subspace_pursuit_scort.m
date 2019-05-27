function [iiter,Lp,trimP,Xtrim] = subspace_pursuit_scort(alg,C,Gn,Gtheta,L,N,patchno,subdiv,stop,t1,t2,u,Y,dub,SNR,currstr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
% Edited by Pavitra Krishnawamy (pavitrak@mit.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
flag = 0;
count = 0;
t = t2-t1+1;

%% Target Sparsity Level if Doubling Heuristic
if subdiv > 1 && subdiv <= 4    %if we are within cortical ico hierarchy
    if dub == 1                 %if doubling heuristic is 1
        L = 2*length(patchno);  %update target sparsity L with "doubling heuristic"
    end
end

%% Compute Adjacency Matrices (Adj) based on Mutual Coherence Thresholds (u)
Gcorr = abs(Gn'*Gn);
Adj = Gcorr > u(subdiv);

%% Greedy Subspace Pursuit Algorithm
while isequal(flag,0)
    count = count + 1;
    Xtrim{count} = zeros(size(Gtheta,2),t);
    Xtrim_us{count} = zeros(size(Gtheta,2),t);
    
    if  isequal(count,1)
        % Compute initial proxy and sort based on L2-norm
        [proxy, proxy_us] = source_estimates(alg,C,Gtheta,N,Y,t1,t2,SNR,stop,currstr);
        [~,IND] = sort(sqrt(dot(proxy,proxy,2)/t),'descend');
        
        % Initialize conflict vector (CV) to refine support (IX)
        l = 1;
        IX{count} = IND(l);
        CV = Adj(:,IX{count});
        count2 = l;
        
        while count2 < L && l < length(IND)
            l = l + 1;
            % Check for mutual coherence or adjacency conflicts
            if isequal(CV(IND(l),1),0)
                % Concatenate support indices and update CV
                IX{count} = cat(2,IX{count},IND(l));
                CV = CV + Adj(:,IND(l));
                count2 = count2 + 1;
            end
        end
        
        % Update L based on size of refined support
        Lp = count2;
        
        % Expectation over trimmed support indices (trimX)
        trimX{count} = IX{count}; disp(trimX{count})
        [Xtrim{count}(trimX{count},:), Xtrim_us{count}(trimX{count},:)] = source_estimates(alg,C,Gtheta(:,trimX{count}),N,Y,t1,t2,SNR,stop,currstr(trimX{count}));
        figure; plot(sqrt(dot(Xtrim{count},Xtrim{count},2)));
        
        % Compute residual (rEs) and residual energy (E)
        rEs{count} = Y - Gtheta*Xtrim_us{count};
        E(count) = norm(rEs{count},'fro')^2;
        
    else
        % Compute proxy on residual and sort to ID undiscovered components
        [proxy, proxy_us] = source_estimates(alg,C,Gtheta,N,rEs{count-1},t1,t2,SNR,stop,currstr);
        [~,IND] = sort(sqrt(dot(proxy,proxy,2)/t),'descend');
        
        % Expand support
        ix = [];
        l = 0;
        count2 = l;
        CV = sum(Adj(:,trimX{count-1}),2);
        
        while count2 < Lp && l < length(IND) % count2 < 2*Lp && l < length(IND)
            l = l + 1;
            % Check for adjacency conflicts
            if isequal(CV(IND(l),1),0)
                % Concatenate support indices and update CV
                ix = cat(2,ix,IND(l));
                CV = CV + Adj(:,IND(l));
                count2 = count2 + 1;
            end
        end
        
        % Expectation over expanded support
        IX{count} = unique(cat(2,trimX{count-1},ix));
        [Xls, Xls_us] = source_estimates(alg,C,Gtheta(:,IX{count}),N,Y,t1,t2,SNR,stop,currstr(IX{count}));
        
        % Sort estimates and trim support
        [~,IND2] = sort(sqrt(dot(Xls,Xls,2)/t),'descend');
        trimX{count} = IX{count}(1,IND2(1:Lp)); disp(trimX{count});
        
        % Expectation over trimmed support
        [Xtrim{count}(trimX{count},:), Xtrim_us{count}(trimX{count},:)] = source_estimates(alg,C,Gtheta(:,trimX{count}),N,Y,t1,t2,SNR,stop,currstr(trimX{count}));
        %figure; plot(sqrt(dot(Xtrim{count},Xtrim{count},2)));
        
        % Compute residual and residual energy
        rEs{count} = Y - Gtheta*Xtrim_us{count};
        E(count) = norm(rEs{count},'fro')^2;
        
        % Check for convergence before moving on to next iteration
        [min_E, e_i] = min(E);
        
        for ii = 1:count-1
            s = setdiff(trimX{count},trimX{ii});
            
            if isempty(s)   % if converged, check iff minimum residual
                flag = 1;
                iiter = length(trimX);
                %save('test.mat','trimX','Xtrim');
                if E(count) <= min_E
                    disp('converged per usual');
                    trimP = trimX{ii};
                    trimX{count} = trimX{ii};
                    Xtrim{count} = Xtrim{ii};
                else
                    disp('converged per residual error min.');
                    trimP = trimX{e_i};
                    trimX{count} = trimX{e_i};
                    Xtrim{count} = Xtrim{e_i};
                end
                %figure; plot(sqrt(dot(Xtrim{count},Xtrim{count},2)));
                break  
            end
        end
    end
end
%figure, set(gcf,'color','white'); plot(E); %pause; 
%close all;
end