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
% ***********************************************************************
function [roi] = mia_stats_permutations_rois(roi,OPTIONS) 

%  Loop through regions
for ii=1:length(roi) 

      %% Get current roi signals 
      croi = roi{ii};
      sig1 = croi.sig1;
      sig2 = croi.sig2;
         
for pp=1:OPTIONS.nperm 
    
    for cc=1:size(sig1,1)
    
        % Collect the trials of the two experimental condition in a single
        % set 
        d = cat(3,sig1{cc},sig2{cc});
        n = randperm(size(d,3)); 
        ff1{cc} = d(:,:,n(1:size(sig1{cc},3))) ; 
        ff2{cc} = d(:,:,n(size(sig1{cc},3):end)) ; 
        
    end

        
    % Compute unpaired ttest 
    [tvals,~]=compute_ttest2distrib(cat(3,ff1{:}),cat(3,ff2{:}));
   
    % Smooth the tvalues
    tvalsf = filtfilt(ones(1,OPTIONS.smoth)/OPTIONS.smoth,1,tvals);
    
    % Compute threshold for tvalues (df = nSample(CondA) + nSample(condB) -2)
    thresht = tinv(1-OPTIONS.threshp/2,size(cat(3,ff1{:}),3)+size(cat(3,ff2{:}),3)-2);
    
    % Get the positive segment on tvalues vector
    [seg_p,r_edge_p] = get_significant_segments(tvalsf>thresht); 
    
    % Get the negative segment on tvalues vector
    [seg_n,r_edge_n] = get_significant_segments(tvalsf<-thresht); 
    
    % Get all segment (positive and negative)
    [seg,r_edge] = get_significant_segments(abs(tvalsf)>thresht); 
        
       % For positive clusters
    if isempty(seg_p) 
       max_sumt(pp) = 0 ; 
    else

        % For all segments compute sum of tvals
        for ss=1:length(seg_p)
            sumt_p(ss) = sum(tvalsf(r_edge_p(ss):r_edge_p(ss)+seg_p(ss)));
        end
        max_sumt(pp) = max(sumt_p);
    end
   
    % For negative clusters
    if isempty(seg_n)
       min_sumt(pp) = 0 ; 
    else 
        % For all segments compute sum of tvals
        for ss=1:length(seg_n)
            sumt_n(ss) = sum(tvalsf(r_edge_n(ss):r_edge_n(ss)+seg_n(ss)));
        end
        min_sumt(pp) = min(sumt_n);
    end
   
   
   % for thresholding on the duration
   if isempty(seg)
        max_dur(pp)  = 0 ; 
    else
        max_dur(pp) = max(seg); 
    end
    
end
    
    % Save the results
    roi{ii}.thresh_dur = quantile(max_dur,0.95);
    roi{ii}.thresh_sumtp = quantile(max_sumt,0.95);
    roi{ii}.thresh_sumtn = -quantile(-min_sumt,0.95);
    
    roi{ii}.dur = max_dur;
    roi{ii}.sumtp = max_sumt;
    roi{ii}.sumtn = min_sumt;
    
end


