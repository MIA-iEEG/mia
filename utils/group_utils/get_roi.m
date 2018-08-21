function [roi] = get_roi(m_table_effect,t, s, smask, all_labels,freqs,opt) 
%
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
% Copyright (C) 2016-2018 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
%
% 2015/5/15 : ASD creation

id_loc2 =5; 
id_lat=4;
starting_freq=6  ;
id_subj = 1 ;
roi = [] ;
ptCor = [] ;
roi_rejected = [];
m_bool_res = [] ; 

n  = 2; % Nombre de contact min par roi

% Create unique ROI name
[un,~,ic] = unique(strcat(m_table_effect(:,id_lat),'.',m_table_effect(:,id_loc2)),'stable');

% Hist des occurences
[N,X]= hist(ic,length(un));

idx=find(N>=n) ; 

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
[unsubj, ~,icsubj]=unique(m_table_effect(:,id_subj)) ; 

% Distinct color for subj
clr = hsv(numel(unsubj));  

% ROI counter
ct2 = 1; 

if opt.signifmode ~=0
    bool_signif = cell2mat(m_table_effect(:,starting_freq+opt.freq)) ; 
else
     bool_signif = ones(size(m_table_effect,1),1) ;
end

% For all region that has more than N contacts 
for jj=1:length(idx)

    % Boolean  for one region
    bool_roi = ismember(strcat(m_table_effect(:,id_lat),'.',m_table_effect(:,id_loc2)),un(idx(jj))) ; 

    % If at leaste two patients have active contacts 
    subj_active =unique(m_table_effect(logical(bool_roi.*bool_signif),id_subj)) ; 
  
    if length(subj_active)>=opt.nPt

        str_roi = cell2mat(un(idx(jj))) ; 
       
        % Initialization
        ct = 1;
        subj_in = [];
        masked_sig = [] ;
        labels_roi = [] ;
        m_res_roi = [] ; 
        m_bool_roi = [] ; 
        all_sig = [] ; 
        idpt = [] ; 
        
        % Loop over all subjs
        for ii=1:length(unsubj)

            % All contacts present a region for a pt
            ind_roi_subj = logical((ii==icsubj).*bool_roi.*bool_signif); 

            % If this suject has an activation in this roi
            if sum(ind_roi_subj~=0)
                
                % Subj idx
                subj_in(ct) = ii ; 
                
                % Filter out double occurences
                [c,ia,ic] =  unique(m_table_effect(ind_roi_subj,6)) ; 
                idx1 = find(ind_roi_subj~=0);
                
                % Pick signals
                sig = s(idx1(ia),:,opt.freq)';
                
                % Pick significativity masks
                masked_sig = cat(2,masked_sig,smask(idx1(ia),:,opt.freq)');
                
                % Pick labels
                current_labels = all_labels(idx1(ia)) ; 
                labels_roi = cat(1,labels_roi,current_labels);
             
                m_bool_roi = cat(1,idx1(ia),m_bool_roi);
                
                % Mean signals per patients
                if strcmp(opt.signmode,'signed') 
                    mean_sig_subj(:,ct) = mean(sig,2) ; 
                else 
                    mean_sig_subj(:,ct) = mean(abs(sig),2) ; 
                end
                
                % Concatenates all contacts signals
                all_sig = cat(2,all_sig,sig);
                idpt = cat(1, idpt, repmat(ct,size(sig,2),1)) ;                 
        
                % Increase subj indice
                ct=ct+1; 
    
            end

        end
        
        % Find onset of the first significant period
        d = sum(masked_sig',1) ; 
        detect = diff(d) ; 
        onset = find(detect~=0,1) ; 
        roi{ct2}.onset= t(onset) ; 
 
        % Compute correlations between contacts
        r = corrcoef(all_sig) ; 
        
        % If FLIP option is used 
        if isfield(opt, 'flip_thresh')
           
           % Select min in R 
           A = (r<opt.flip_thresh) ; 
           Ap =  (r>-opt.flip_thresh) ; 
           
           % Finds ind with anti-correlation 
           [Y,I] = max(sum(A)) ;              
           fl = ones(1,length(A)) ; 
           fl(find(Ap(I,:))) = -1 ; 
             
           %flip other signs to positively
           % correlate with the frst timserie
           all_sigFl= all_sig.*fl;    
             
           % Recompute correlations between FLIPPED contacts
           r = corrcoef(all_sigFl) ; 
           all_sig = all_sigFl ;
           
           % Recompute means poer patients
           for ss=1:max(idpt) 
                mean_sig_subj(:,ss) = mean(all_sig(:,idpt==ss),2);
           end
           
           masked_sig = masked_sig.*fl ;
        
           % Add _FLP at the end of the contacts labels that were flipped
           labels_roi(fl==-1) = strcat(labels_roi(fl==-1),'_FLP') ; 

        end
        
        % Compute correlations between patients : interpatient correlation
       rPt = corrcoef(mean_sig_subj) ; 
       trPt = tril(rPt,-1) ; 
       roi{ct2}.corrPt = mean(trPt(trPt~=0)); 
       roi{ct2}.allCorrPt = trPt(trPt~=0);        
       
       % Compute correlations between channels : intra+inter patients
       trChan = tril(r,-1) ; 
       roi{ct2}.corrChan =  mean(trChan(trChan~=0)); 
       roi{ct2}.allCorrPt = trChan(trChan~=0);
       
       roi{ct2}.name =  cell2mat(un(idx(jj))) ; 
       roi{ct2}.signmoy = mean_sig_subj ;
       roi{ct2}.Fmask = masked_sig ;
       roi{ct2}.labels = labels_roi ;
       roi{ct2}.idPt = subj_in ; 
       roi{ct2}.namePt = unsubj(subj_in) ; 
       roi{ct2}.t = t ; 
       roi{ct2}.freq= freqs{opt.freq} ; 
       
       ct2=ct2+1;

        clear rC trC ptCor chanCor r tr mean_sig_subj A ; 

 
    end

end


