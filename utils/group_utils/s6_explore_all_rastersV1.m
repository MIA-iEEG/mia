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
function [roi] = s6_explore_all_rastersV1(roi,datafiles,RTfiles,OPTIONS) 

ct=1;
%  Loop through regions
for ii=1:length(roi) 

    croi = roi{ii};
    
    k=strfind(croi.labels,'_');
    
    for pp=1:length(k)
    
        str = cell2mat(croi.labels(pp)) ;
        
        % Column for pt ID
        if strcmpi(OPTIONS.mtg,'monopolar') ; tt = 0 ; else tt=1 ; end
        
        m_tab{ct,1} = str(1:k{pp}(end-tt)-1);
        
        % Column for electrode
        m_tab{ct,2} = str(k{pp}(end-tt)+1:end);
        
        % Add column for indexing the roi
        m_tab{ct,3} = ii;
        ct=ct+1;
        
    end
end

% Gets all unique occurence of patients
[un, ia, ic] = unique(m_tab(:,1),'stable');


%reorder data files
for ii=1:length(un)
    idx = ~logical(cellfun(@isempty,strfind(datafiles,un{ii})));
    dfiles(ii,1)=datafiles(idx);
end
datafiles=dfiles;
clearvars dfiles

%reorder RT files
for ii=1:length(un)
    idx = ~logical(cellfun(@isempty,strfind(RTfiles,un{ii})));
    rtFILES(ii,1)=RTfiles(idx);
end
RTfiles=rtFILES;
clearvars rtFILES
                
rasters = cell(size(m_tab,1),1);
rts = cell(size(m_tab,1),1);

for kk=1:length(un)
     
if isempty(strfind(RTfiles{kk},'MiSsInG'))

    dat = load(datafiles{kk});      

    dat_rt=load(RTfiles{kk});
    
    
     
    isGood = (dat_rt(find(dat_rt(:,4)),3)==1)&dat_rt(find(dat_rt(:,4)),2)<(OPTIONS.win_noedges(2)*1000);
    
    isGoodinRt = (dat_rt(:,3)==1)&(dat_rt(:,2)<(OPTIONS.win_noedges(2)*1000))&dat_rt(:,4);

    % Get channel to display
    idx_elec_tab = ismember(m_tab(:,1),un{kk}) ;
    lab_toproc = m_tab(idx_elec_tab,2);
    idx_elec_dat = ismember(dat.labels,lab_toproc);    

    
    id = find(idx_elec_dat);
  
    [LIA,LOCB] = ismember(lab_toproc',dat.labels(idx_elec_dat)) ;
    
    % Format output
    rasters(idx_elec_tab) =  num2cell(dat.zs(id(LOCB),(dat.Time>OPTIONS.win_noedges(1))&(dat.Time<OPTIONS.win_noedges(2)),isGood),[2,3]);
    rts(idx_elec_tab) = num2cell(dat_rt(isGoodinRt,2),1) ;
    op(idx_elec_tab) = num2cell(dat_rt(isGoodinRt,1),1) ;
    idpt(idx_elec_tab) = repmat({repmat(kk,sum(isGoodinRt),1)},sum(idx_elec_tab),1);
      
end
end

%     rasters=rasters(~logical(cellfun(@isempty,resters)));
%     rts=rts(~logical(cellfun(@isempty,rts)));
%     op=op(~logical(cellfun(@isempty,op)));
%     idpt=idpt(~logical(cellfun(@isempty,idpt)));
    
    
% Reinject the signals into the roi structure 
for ii=1:length(roi) 

    idx = [m_tab{:,3}] == ii ; 
    
    roi{ii}.F =rasters(idx);
    roi{ii}.rts = rts(idx);
    roi{ii}.op= op(idx);
    roi{ii}.idptsignals= idpt(idx);
    roi{ii}.missingPt=[];
    
end

%remove patient if TR file is missing

for kk=1:length(un)
    if ~isempty(strfind(RTfiles{kk},'MiSsInG'))
        name2vir=RTfiles{kk};
        %get back missing subject name
        name2vir=name2vir(1:length(strtok(name2vir,'M'))-1);      
        
        %look for the existence of this patient in Roi
        for ii=1:length(roi)
            
            %get indices to remove corresponding datas
            id2vir1=find(strcmp(roi{ii}.namePt,name2vir));            
            id2vir2=find(~cellfun(@isempty,strfind(roi{ii}.labels,name2vir)));
           
            %remove corresponding datas
            if ~isempty(id2vir1) && ~isempty(id2vir2)
                roi{ii}.signmoy(:,id2vir1)=[];
                roi{ii}.idPt(id2vir1)=[];
                roi{ii}.namePt(id2vir1)=[];
                
                roi{ii}.Fmask(:,id2vir2)=[];
                roi{ii}.labels(id2vir2)=[]; 
                roi{ii}.F(id2vir2)=[]; 
                roi{ii}.rts(id2vir2)=[]; 
                roi{ii}.op(id2vir2)=[]; 
                roi{ii}.idptsignals(id2vir2)=[]; 
                
                %create a field to keep the information of missing patient
                %(could be display in figure ??)
                if ~isempty(roi{ii}.missingPt)                    
                    roi{ii}.missingPt{length(roi{ii}.missingPt)+1}=name2vir;
                else
                    roi{ii}.missingPt{1}=name2vir;
                end
            end
            
            
    
        end
    end
end




























