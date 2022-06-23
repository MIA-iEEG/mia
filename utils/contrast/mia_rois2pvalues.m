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
function [lroi] = mia_rois2pvalues(roi,OPTIONS) 

%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    
    %% Get current roi signals 
    sig1 = croi.sig1;
    sig2 = croi.sig2;
    mat_sig1 = squeeze(cat(3,sig1{:})) ; 
    mat_sig2 = squeeze(cat(3,sig2{:})) ; 
    sig1cat = cat(3,sig1{:});
    sig2cat = cat(3,sig2{:});
   
    % Compute TTEST
    [tvals,~]=compute_ttest2distrib(sig1cat,sig2cat);
    
    tvalsf = filtfilt(ones(1,OPTIONS.smoth)/OPTIONS.smoth,1,tvals);
    thresht = tinv(1-OPTIONS.threshp/2,size(sig1cat,3)+size(sig2cat,3)-2);    
    
    % Duration
    [seg,r_edge] = get_significant_segments(abs(tvalsf)>thresht); 
  
    [Y,I] = sort(croi.dur);
    pv = [] ;
    Fs = 1 /(croi.t(2) - croi.t(1));
    
    for ss=1:length(seg) 
        pv(ss) = length(find(Y>seg(ss)))/Fs;
    end
    
    stats.pdur{ii}= pv ;
    
    %% Get clusters of positive and negative t-values (separately)
    [seg_p,r_edge_p] = get_significant_segments(tvalsf>thresht); 
    [seg_n,r_edge_n] = get_significant_segments(tvalsf<-thresht); 
    [Ypos,I] = sort(croi.sumtp);
    [Yneg,I] = sort(croi.sumtn);
    
    %% Process positive clusters
    pvpos = [] ;

    % For all segments (positive) compute sum of tvals
    for ss=1:length(seg_p)
        sumt_p = sum(tvalsf(r_edge_p(ss):r_edge_p(ss)+seg_p(ss)));
        pvpos(ss) = length(find(Ypos>sumt_p))/Fs;

    end
    stats.psumtp{ii}=pvpos; 
    
    %% Process negative clusters
    pvneg = [] ;
     % For all segments (negative) compute sum of tvals
     for ss=1:length(seg_n)
        sumt_n = sum(tvalsf(r_edge_n(ss):r_edge_n(ss)+seg_n(ss)));
        pvneg(ss) = length(find(Yneg<sumt_n))/Fs;
     
     end
     stats.psumtn{ii}=pvneg ; 
    
     % Create a result structure 
     lroi{ii}.name = croi.name; 
     lroi{ii}.signmoy = croi.signmoy; 
     lroi{ii}.idPt = croi.idPt ;
     lroi{ii}.tvalsf= tvalsf ;
     lroi{ii}.thresht= thresht ;
     lroi{ii}.meansig1= mean(mat_sig1,2) ; 
     lroi{ii}.meansig2= mean(mat_sig2,2) ; 
     lroi{ii}.t = croi.t;
     lroi{ii}.psumtn = pvneg;
     lroi{ii}.psumtp = pvpos ; 
     lroi{ii}.pdur = pv ; 
     lroi{ii}.namePt = croi.namePt;
     lroi{ii}.meancorr = croi.meancorr;
     lroi{ii}.onsetdur = croi.t(r_edge);
     lroi{ii}.onsetpvalp = croi.t(r_edge_p);
     lroi{ii}.onsetpvaln = croi.t(r_edge_n);

end


