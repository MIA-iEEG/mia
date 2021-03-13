function [roi] = get_roi(m_table_effect,t, s, smask, all_labels,freqs,opt) 
%
% ***********************************************************************
%
%  Copyright (C) 2016-2021 CNRS - Universite Aix-Marseille
%
%  This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
%
% ***********************************************************************
%
% 2015/5/15 : ASD creation
% 2017 : ASD modif flip according to description : 
% We computed the correlation matrix across all signals within one region 
%   (i.e. across contacts and patients).
%   This matrix was thresholded separately for positive and negative values, 
%   based on a threshold chosen by the user (default thr = -0.7). 
%   We counted the number of high negative correlation (-thr) 
%   and pick the channel with the highest number of anti-correlated channels.
% All channels highly correlated with this channel were flipped. 
%
% 2018/12/7 : correct flip : do not flip if nothing needs to be flipped! 
%
% 2019/10/28 : ASD : remove idx to count region - useless

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

% Show all contacts or Show only significant
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
 
    all_sig = s(idx_signals,:,opt.freq)';
    masked_sig = smask(idx_signals,:,opt.freq)';
    labels_roi = all_labels(idx_signals);
    subj_in = find(ismember(unsubj,subj_active)); 
 
    % Get each contact's patient name 
     for cc=1:length(labels_roi)
        ptchar = labels_roi{cc} ; 
        idx_underscores = strfind(ptchar,'_');
        %ptname{cc} = ptchar(1:idx_underscores(2)-1);
        %ptname{cc} = ptchar(1:idx_underscores(1)-1);
        if strcmp(opt.montage,'bipolar')
            ptname{cc} = ptchar(1:idx_underscores(end-1)-1);
        else 
            ptname{cc} = ptchar(1:idx_underscores(end)-1);
        end
     end    
    
%     idx_underscores = cellfun(@(x) strfind(x,'_'),labels_roi, 'UniformOutput',false) ; 
%     idx_patient = cell2mat(cellfun(@(x) x(end-1), idx_underscores,'UniformOutput',false)) ; 
    
 % If FLIP option is used 
    if isfield(opt, 'flip_thresh')
          [all_sig,masked_sig,labels_roi,r] = flip_signals(all_sig, masked_sig, labels_roi, opt.flip_thresh) ; 
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
    trPt = tril(rPt,-1) ; 
    roi{ctRoi}.corrPt = mean(trPt(trPt~=0)); 
    roi{ctRoi}.allCorrPt = trPt(trPt~=0);        

    % Compute correlations between channels : intra+inter patients       
    trChan = tril(r,-1) ; 
    roi{ctRoi}.corrChan =  mean(trChan(trChan~=0)); 
    roi{ctRoi}.allCorrChan = trChan(trChan~=0);
   
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
function  [all_sig,masked_sig,labels_roi,r]  =  flip_signals(all_sig, masked_sig, labels_roi, flip_thresh)

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

%     % Recompute means per patients
%     [C,IA,IC] = unique(cellfun( @(x) x(1:3), labels_roi, 'UniformOutput',false )) ;
%     for ss=1:max(IC) 
%         mean_sig_subj(:,ss) = mean(all_sig(:,IC==ss),2);
%     end
    
    % Add _FLP at the end of the contacts labels that were flipped
    labels_roi(fl==-1) = strcat(labels_roi(fl==-1),'_FLP') ; 
