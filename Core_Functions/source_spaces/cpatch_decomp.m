function cpatch_decomp(subjid, fwd, whiteG, nthresh, p_sing, pt_s, subdiv_stop, svds_fname, meastype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE SINGULAR VALUE DECOMPOSITION (SVD) OF DISJOINT PATCH LEAD FIELDS
% COMPUTE NORMALIZED MEAN REPRESENTATION ACCURACY (NMRA) TO PICK EFFECTIVE NUMBER OF MODES PER PATCH
% Written by Gabriel Obregon-Henao (obregon@nmr.mgh.harvard.edu)
% Modified by Pavitra Krishnaswamy (pavitrak@nmr.mgh.harvard.edu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute total number of dipoles in left hemisphere of densely-sampled source space
nuse = double(fwd.src(1).nuse);

%% Load vertex numbers of superficial dipoles (< 5mm to inner_skull.surf) 
%Left hemisphere
fid = fopen(strcat(subjid,'-meg-all-mindistout-lh-label.txt'));
cols = textscan(fid,'%f');
fclose(fid);
lk_s = cols{1};

%Right hemisphere 
fid = fopen(strcat(subjid,'-meg-all-mindistout-rh-label.txt'));
cols2 = textscan(fid,'%f');
fclose(fid);
lk2_s = cols2{1};

for subdiv = 1:min(subdiv_stop,4)
    %% Initialize
    PT = 0;
    %Load source spaces with cortical patch statistics
    source(subdiv).hemisph = mne_read_source_spaces(strcat(subjid,'-ico-',num2str(subdiv),'p-src.fif'));
    %Total number of patches in left hemisphere
    puse = length(source(subdiv).hemisph(1).pinfo);
    disp(['Computing SVD of disjoint patch lead fields in ico-' num2str(subdiv) 'p source space... '])

    for ii = 1:puse
        %% Form disjoint lead field of patch ii in left hemisphere (Gx)
        [~,~,ico(subdiv).patch(ii).lk] = intersect(source(subdiv).hemisph(1).pinfo{ii},fwd.src(1).vertno,'legacy');
        Gx = whiteG(:,ico(subdiv).patch(ii).lk);
        
        %% Compute SVD of Gx and Keep <p_sing>-largest modes of Gx
        [ico_full.patch(ii).U,ico_full.patch(ii).S,ico_full.patch(ii).V] = svd(full(Gx));
        ico(subdiv).patch(ii).G = Gx;
        %p_sing_use  = min(p_sing, size(ico_full.patch(ii).S,1));
        ico(subdiv).patch(ii).U = ico_full.patch(ii).U(:,1:p_sing);
        ico(subdiv).patch(ii).S = ico_full.patch(ii).S(1:p_sing,1:p_sing);
        ico(subdiv).patch(ii).V = ico_full.patch(ii).V(:,1:p_sing);

        %% Form disjoint lead field of patch ii+puse in right hemisphere (Gx2)
        [~,~,lk2] = intersect(source(subdiv).hemisph(2).pinfo{ii},fwd.src(2).vertno,'legacy');
        ico(subdiv).patch(ii+puse).lk = lk2 + nuse;
        Gx2 = whiteG(:,ico(subdiv).patch(ii+puse).lk);
        
        %% Compute SVD of Gx2 and Keep <p_sing>-largest modes of Gx2
        [ico_full.patch(ii+puse).U,ico_full.patch(ii+puse).S,ico_full.patch(ii+puse).V] = svd(full(Gx2));
        ico(subdiv).patch(ii+puse).G = Gx2;
        %p_sing_use2  = min(p_sing, size(ico_full.patch(ii+puse).S,1));
        ico(subdiv).patch(ii+puse).U = ico_full.patch(ii+puse).U(:,1:p_sing);
        ico(subdiv).patch(ii+puse).S = ico_full.patch(ii+puse).S(1:p_sing,1:p_sing);
        ico(subdiv).patch(ii+puse).V = ico_full.patch(ii+puse).V(:,1:p_sing);
        
        %% Check if patches ii and ii+puse are superficial (at least <pt_s> percent of dipoles are < 5mm to inner_skull.surf)
        [int,~,~] = intersect(ico(subdiv).patch(ii).lk,lk_s,'legacy');
        [int2,~,~] = intersect(lk2,lk2_s,'legacy');
        if round(length(int)/length(ico(subdiv).patch(ii).lk)*100) >= pt_s
            disp(['Patch ' num2str(ii) ' labelled as superficial (at least ' num2str(pt_s) '% of dipoles are < 5mm to inner_skull.surf)'])
            ico(subdiv).patch(ii).mindistout = true(1);
        else
            ico(subdiv).patch(ii).mindistout = false(1);
        end    
        if round(length(int2)/length(lk2)*100) >= pt_s
            disp(['Patch ' num2str(ii+puse) ' labelled as superficial (at least ' num2str(pt_s) '% of dipoles are < 5mm to inner_skull.surf)'])
            ico(subdiv).patch(ii+puse).mindistout = true(1);
        else
            ico(subdiv).patch(ii+puse).mindistout = false(1);
        end
        
        %% Define neighborhoods (nearest neighbors) across entire patch source space
        ico = patch_neighbors(ico,ii,puse,source,subdiv);
        
        %% SVDs computation: percentage completed
        pt = round(100*ii/puse);       
        if ~mod(pt,20) && ~isequal(PT,pt)
            disp(['Computing... ' num2str(pt) ' percent [done]'])
            PT = pt; 
        end
        clear Gx Gx2
    end
    
    %% Find superficial patch numbers (at least <pt_s> percent of dipoles are < 5mm to inner_skull.surf)
    B = [ico(subdiv).patch.mindistout];
    ii_s = find(B);
    clear B
    
    %% NMRA computations in non-superficial patches using <p_sing> modes
    disp(['Computing NMRA using p = 1:',num2str(p_sing),' modes...'])
    figure, set(gcf, 'color', 'white') 
    for ii = 1:p_sing
        % Compute NMRA
        nmra = [];
        for jj = 1:length(ico(subdiv).patch)
            if ~ismember(jj,ii_s)
                diagS = diag(ico_full.patch(jj).S);
                nmra = cat(2,nmra,dot(diagS(1:ii)',diagS(1:ii)',2)/dot(diagS(:)',diagS(:)',2));          
            else
                if isequal(ii,1)
                    disp(['Excluding superficial patch ' num2str(jj) ' from NMRA computations'])
                end
            end
        end
        
        % Plot NMRA Histogram
        subplot(2,round(p_sing/2),ii); hist(nmra)
        grid on; title(['ico-',num2str(subdiv),'; p = ',num2str(ii),' modes'])
        xlabel('NMRA'); ylabel('Number of patches')
        
        % Plot Median Energy captured by <p_i> Modes
        E(ii) = median(nmra);
        hold on; plot(E(ii),0,'*r')
    end
    h = gcf; print(h,'-dpdf',strcat(subjid, '_nmras_',meastype,'_ico',num2str(subdiv),'_prewhite'));
    
    %% Choose effective number of modes per patch
    B2 = find(round(E*100) >= nthresh);
    p(subdiv) = B2(1);
    clear ico_full
    
end

%% Save File
save(svds_fname,'ico','source','p','-v7.3')
clear ico source 

%% Add on Patch Areas to SVDs File
patch_areas(subjid, svds_fname, fwd, subdiv_stop); 
h = gcf; print(h,'-dpdf',strcat(subjid, '_pareas_',meastype,'_prewhite'));
    
end