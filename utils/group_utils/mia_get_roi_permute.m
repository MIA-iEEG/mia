% ========================================================================
% This file is part of MIA.
% 
% MIA is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MIA is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%  
% Copyright (C) 2016-2022 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
function [roi] = mia_get_roi_permute(m_table_effect,t, s, smask, all_labels,freqs,opt) 

id_loc2 =5; 
id_lat=4;
starting_freq=6  ;
id_subj = 1 ;
roi = [] ;

% Create unique ROI name
[un,~,ic] = unique(strcat(m_table_effect(:,id_lat),'.',m_table_effect(:,id_loc2)),'stable');

if ~isempty(opt.freq)
    
    % Loop could be remove, for display only..
    for jj=1:length(un)

        idx_tmp = (ic==jj)&(cell2mat(m_table_effect(:,starting_freq+opt.freq))==1);
        idx_tmp2 = (ic==jj);
        
        N_all = length(unique(m_table_effect(idx_tmp2,6)));
        
        % Get unique occurence of electrodes in m_table_effect
        N_active = length(unique(m_table_effect(idx_tmp,6)));

        N_pt = length(unique(m_table_effect(idx_tmp,1))) ;

        fprintf(sprintf('%s \t %d/%d \t contacts \t(%d pt) \n',un{jj},N_active,N_all,N_pt));

    end
end

% Get unique occurence of patients
[unsubj, ~,~]=unique(m_table_effect(:,id_subj)) ; 

% Distinct color for subj
clr = hsv(numel(unsubj));  

% ROI counter
ctRoi = 1; 

% Filters contacts (or not) based on significance of signal
if opt.signifmode ~=0
    bool_signif = cell2mat(m_table_effect(:,starting_freq+opt.freq)) ; 
else
    bool_signif = ones(size(m_table_effect,1),1) ;
end

% For all unique regions 
for jj=1:length(un)

    % Boolean for one region
    bool_roi = ismember(strcat(m_table_effect(:,id_lat),'.',m_table_effect(:,id_loc2)),un(jj)) ; 

    % Gets patients which contribute to this region
    subj_active =unique(m_table_effect(logical(bool_roi.*bool_signif),id_subj)) ; 
   
    % Do not take into account regions which have less patients than nPt 
    if length(subj_active)<opt.nPt ; continue ; end
   
    % Removes doublons and gets contacts : involved in ROI / significant
    [c,ia,ic] =  unique(m_table_effect(logical(bool_roi.*bool_signif),6)) ; 
    
    % Gets indices with doulbons
    idx_signals_doublons = find(logical(bool_roi.*bool_signif)) ; 
    
    % Gets indices without doublons
    idx_signals = idx_signals_doublons(ia) ; 
 
    % PERMUTATIONS
    for pp = 1:1000

        idx_signals_perm = [] ; 

        % Loop over patients with active contacts in this ROI
        for ii=1:length(subj_active)
            
            % Current patient boolean in m_table_effect
            bool_current_pt  = ismember(m_table_effect(:,id_subj),subj_active(ii)) ; 

            % All contacts present for this region and this pt
            ind_roi_subj = logical(bool_current_pt.*bool_roi.*bool_signif); 

            % Here for permutation we select random contacts from this patient 
            nbToSelect = length(unique(m_table_effect(ind_roi_subj,6)));
            [c2,ia2,ic2] =  unique(m_table_effect(bool_current_pt,6)) ; 
            
            n = randperm(length(ia2),nbToSelect);

            idx_signals_perm = cat(2,idx_signals_perm,ia2(n)'+find(bool_current_pt,1)-1) ; 

            % Get indices of signals for this region/patient (valid in
            % m_table_effect, s, smask and all_labels)
    %         idx_signals = cat(1,idx_signals,idx1(ia)) ;      

        end

        all_sig = s(idx_signals_perm,:,opt.freq)';

        % If FLIP option is used 
        if isfield(opt, 'flip_thresh')
            [~,~,labels_roi,r] = flip_signals(all_sig, smask(idx_signals_perm,:,opt.freq)', all_labels(idx_signals_perm), opt.flip_thresh) ; 
        else
            r = corrcoef(all_sig) ;
            labels_roi = all_labels(idx_signals_perm) ; 
        end

        % Compute correlations between contacts
        r = atanh(r);                                   % apply Fisher-tansformation (convert r to Fisher's z)
        trChan = tril(r,-1) ; 
        corr_permut(pp) = tanh(mean(trChan(trChan~=0))); % Compute the average of all Fisher's z values from the lower triangular part then convert to r value 
        labels_permut{pp} = labels_roi;
        
        if length(unique(all_labels(idx_signals_perm))) ~= length(all_labels(idx_signals_perm))
            fprintf('Doublons exist\n') ;
        end

    end


    all_sig = s(idx_signals,:,opt.freq)';
    masked_sig = smask(idx_signals,:,opt.freq)';
    labels_roi = all_labels(idx_signals);
    subj_in = find(ismember(unsubj,subj_active)); 
 
    
        % Get each contact's patient name 
     for cc=1:length(labels_roi)
        ptchar = labels_roi{cc} ; 
        idx_underscores = strfind(ptchar,'_');
        if strcmp(opt.montage,'bipolar')
            ptname{cc} = ptchar(1:idx_underscores(end-1)-1);
        else 
            ptname{cc} = ptchar(1:idx_underscores(end)-1);
        end
     end  

     % If FLIP option is used 
    if isfield(opt, 'flip_thresh')
          [mean_sig_subj,masked_sig,labels_roi,r] = flip_signals(all_sig, masked_sig, labels_roi, opt.flip_thresh) ; 
    else
           r = corrcoef(all_sig) ; 
 
    end

    % Compute means per patients
    [~,~,IC] = unique(ptname);
        
    for ss=1:max(IC) 
        % Mean signals per patients
        if strcmp(opt.signmode,'signed') 
            mean_sig_subj(:,ss) = mean(all_sig(:,IC==ss),2); 
        else 
            mean_sig_subj(:,ss) = mean(abs(all_sig(:,IC==ss)),2) ; 
        end
     end

   % Compute correlations between patients : interpatient correlation
   rPt = corrcoef(mean_sig_subj) ; 
   rPt = atanh(rPt);                       % apply Fisher-tansformation (convert r to Fisher's z)
   trPt = tril(rPt,-1) ; 
   roi{ctRoi}.corrPt = tanh(mean(trPt(trPt~=0))); % Compute the average of all Fisher's z values from the lower triangular part then convert to r value
   roi{ctRoi}.allCorrPt = trPt(trPt~=0);        

   % Compute correlations between channels : intra+inter patients   
   r = atanh(r);                                   % apply Fisher-tansformation (convert r to Fisher's z)
   trChan = tril(r,-1) ; 
   roi{ctRoi}.corrChan =  tanh(mean(trChan(trChan~=0))); % Compute the average of all Fisher's z values from the lower triangular part then convert to r value
   roi{ctRoi}.allCorrChan = trChan(trChan~=0);   
   
   % Saves permutations
    roi{ctRoi}.corr_permut = corr_permut ; 
    roi{ctRoi}.labels_permut = labels_permut ; 

    % Find onset of the first significant period
    d = sum(masked_sig',1) ; 
    detect = diff(d) ; 
    onset = find(detect~=0,1) ; 
    roi{ctRoi}.onset= t(onset) ; 

    % All other fields to roi 
    roi{ctRoi}.name =  cell2mat(un(jj)) ; 
    roi{ctRoi}.signmoy = mean_sig_subj ;
    roi{ctRoi}.Fmask = masked_sig ;
    roi{ctRoi}.labels = labels_roi ;
    roi{ctRoi}.idPt = subj_in ; 
    roi{ctRoi}.namePt = unsubj(subj_in) ; 
    roi{ctRoi}.t = t ; 
    roi{ctRoi}.freq= freqs{opt.freq} ; 

    ctRoi=ctRoi+1;

    clear r mean_sig_subj ptname ; 

 
end

% This function flip signals if option is checked (for broadband "LFP"signals) 
function  [mean_sig_subj,masked_sig,labels_roi,r]  =  flip_signals(all_sig, masked_sig, labels_roi, flip_thresh)

    % Computes correlations
    r = corrcoef(all_sig) ; 

    % Select min in R 
    A = (r<flip_thresh) ; 
    Ap =  (r>-flip_thresh) ; 

    % Finds ind with anti-correlation 
    [Y,I] = max(sum(A)) ;

    % Nothing to flip
    if Y ==0; I = []; end

     % Create a vector which indicate if flip is done (-1) or not (+1) 
    fl = ones(1,length(A)) ; 
    fl(find(Ap(I,:))) = -1 ; 

    % Flip signals
    all_sigFl= all_sig.*fl;    
    
    % Flip masks of significance
    masked_sig = masked_sig.*fl ;

    % Recompute correlations between FLIPPED contacts
    r = corrcoef(all_sigFl) ; 
    all_sig = all_sigFl ;

    % Recompute means per patients
    [C,IA,IC] = unique(cellfun( @(x) x(1:5), labels_roi, 'UniformOutput',false )) ;% BUG To fix : 1:5 is length of pt_name
    for ss=1:max(IC) 
        mean_sig_subj(:,ss) = mean(all_sig(:,IC==ss),2);
    end
    
    % Add _FLP at the end of the contacts labels that were flipped
    labels_roi(fl==-1) = strcat(labels_roi(fl==-1),'_FLP') ; 
