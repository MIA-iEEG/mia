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
% 2015/5/15 : ASD creation
function [] = mia_display_rois_conditions_fdr(roi,OPTIONS) 
dur_name ={};
tp_name = {};
tn_name = {};

allroi = [roi{:}];

% % index of duration post-baseline
onsetdur = [allroi.onsetdur] ; 
onsetpval = [allroi.onsetpvalp,allroi.onsetpvaln] ; 

% create boolean vector to identify if cluster are in the tw
twdur = onsetdur>OPTIONS.tw(1)&onsetdur<OPTIONS.tw(2); 
twpval = onsetpval>OPTIONS.tw(1)&onsetpval<OPTIONS.tw(2); 

% FDR : durations
% Compute FDR only for pvalues in tw 
pdur = [allroi.pdur];
indur = pdur(twdur) ;
seg_dur = mia_ft_fdr(indur,OPTIONS.qFDR);

% Get the boolean vector of clusters that were on/off back in the whole tw
twdur(find(twdur)) = seg_dur ; seg_dur = twdur ; 

% FDR : sum(tvlas) POSITIVE VALUES
psumt = [allroi.psumtp,allroi.psumtn] ; 
inp= psumt(twpval); 
seg_sumt = ft_fdr(inp,OPTIONS.qFDR);
seg_sumtp = seg_sumt(1:length([allroi.onsetpvalp]));
seg_sumtn = seg_sumt(length([allroi.onsetpvalp])+1:end);

% Get the boolean vector of clusters that were on/off back in the whole tw
twpval(find(twpval)) = seg_sumt ; seg_sumt = twpval ; 

ctdur = 1;
ctsumtp = 1 ; 
ctsumtn = 1; 


fprintf('\nDURATION : %d input FDR : %d output FDR',length(indur),sum(seg_dur));
fprintf('\nSUMTP: %d input FDR : %d output FDR',length(inp),sum(seg_sumt));


%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    
    % Duration
    [seg,r_edge] = get_significant_segments(abs(croi.tvalsf)>croi.thresht); 
    hdur = zeros(1,length(croi.tvalsf));

    for ss=1:length(seg)
        if seg_dur(ctdur)==1
            hdur(r_edge(ss):r_edge(ss)+seg(ss)) = 1;
            dur_name = cat(1,croi.name,dur_name) ; 
        end
        ctdur = ctdur+1;
    end
    
    % Sum Tvals
    [seg_p,r_edge_p] = get_significant_segments(croi.tvalsf>croi.thresht);
    [seg_n,r_edge_n] = get_significant_segments(croi.tvalsf<-croi.thresht); 

    hsumt = zeros(1,length(croi.tvalsf));
     
    % For all segments compute sum of tvals
     for ss=1:length(seg_p)
        if seg_sumtp(ctsumtp)==1
            hsumt(r_edge_p(ss):r_edge_p(ss)+seg_p(ss)) = 1;
            tp_name = cat(1,croi.name,tp_name) ; 
        end
        ctsumtp = ctsumtp+1;
     end
        
    % For all segments compute sum of tvals
     for ss=1:length(seg_n)
        if seg_sumtn(ctsumtn)==1
            hsumt(r_edge_n(ss):r_edge_n(ss)+seg_n(ss)) = 1;
            tn_name = cat(1,croi.name,tn_name) ; 
        end
        ctsumtn = ctsumtn+1;
     end
    
     if isfield(OPTIONS,'idxToDisp')
         if sum(ismember(OPTIONS.idxToDisp,ii))~=0
             
            if strcmp(croi.name(end),'L')
                % Displays one figure per ROI 
                hfig(ii) = figure('Name',croi.name,'Units','Normalized','Position',[0 , 0, 0.5,0.9]); 
            else 
                % Displays one figure per ROI 
                hfig(ii) = figure('Name',croi.name,'Units','Normalized','Position',[0.5 , 0, 0.5,0.9]); 
    end
    
    % window background white 
    set(gcf,'color','white');
    
    %% Plot mean signals
    ah1 = subplot(5,1,1) ;   
%     set(gca, 'ColorOrder',OPTIONS.clr(croi.idPt,:))
    
    for jj=1:length(croi.idPt)  
        % Plot patient averaged timeseries
        plot(croi.t,croi.signmoy(:,jj)','color',OPTIONS.clr(croi.idPt(jj),:), 'LineWidth',1.03); hold on ;  
    end    
    
    % Add patient' id in the legend (if they are present in the plot!)
    hleg = legend(strrep(croi.namePt,'_','-'),'Location','NorthWest');
    pos = get(ah1,'position'); 
    posleg = get(hleg,'position'); 
    set(hleg,'position',[pos(1)-0.12 pos(2) posleg(3:4)]);
    title(sprintf('%s R = %f',croi.name,croi.meancorr),'FontSize', 30); grid on ;
    grid on ; xlim(OPTIONS.win_noedges);   

    %% Display the two conditions means (accross patients)
    ah3 = subplot(5,1,2) ; 
  
    plot(croi.t,croi.meansig1,'LineWidth',1.5, 'color',[0 0 0]); hold on ; plot(croi.t,croi.meansig2,'LineWidth',0.1,  'color',[0 0 0]); 
%     plot(croi.t,croi.meansig1,'LineStyle','-', 'color',[0 0 0],'LineWidth',2); hold on ; plot(croi.t,croi.meansig2,'LineStyle','--',  'color',[0 0 0],'LineWidth',2); 
    grid on ;  xlim(OPTIONS.win_noedges);   
    hleg = legend('HOMOG','HETEROG');
    pos = get(ah3,'position'); 
    posleg = get(hleg,'position'); 
    set(hleg,'position',[pos(1)-0.12 pos(2) posleg(3:4)]);
    
    %% Display tvals uncorrected
    ah4 = subplot(5,1,3) ; 
%     set(0,'DefaultAxesColorOrder',jet(size(sig1,1)));

    plot(croi.t,croi.tvalsf.*(abs(croi.tvalsf)>croi.thresht),'LineWidth',2,'color','b');  ylim([-7 7]);
    grid on ;  xlim(OPTIONS.win_noedges);  title(strcat('Uncorrected p<',num2str(OPTIONS.threshp)));
  
    %% Display tvals corrected (sumt or dur)
    ah5 = subplot(5,1,4) ; 
   
     plot(croi.t,croi.tvalsf.*hdur,'LineWidth',2,'color','b');  ylim([-7 7]);
     grid on ;  xlim(OPTIONS.win_noedges);   title('Threshold on duration');
     
  
     %% Display tvals corrected (sumt or dur)
    ah6 = subplot(5,1,5) ; 
    
    plot(croi.t,croi.tvalsf.*hsumt,'LineWidth',2,'color','b');  ylim([-7 7]);
        grid on ;  xlim(OPTIONS.win_noedges); title('Threshold on sumt');
        print(sprintf('roi_%s',croi.name),'-dsvg');
    
         end
         
     end
     
end

uniquedur = unique(dur_name) ; 
fprintf('\nDuration :');
for pp=1:length(uniquedur)
    fprintf(strcat(uniquedur{pp},'\t'));
end

uniquet = unique(cat(1,tp_name,tn_name)) ;
fprintf('\nSumtvals  :');
for pp=1:length(uniquet)
    fprintf(strcat(uniquet{pp},'\t'));
end
fprintf('\n')



