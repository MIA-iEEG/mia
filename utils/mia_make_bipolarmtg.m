function [b, bilabels, idx1, idx2] = mia_make_bipolarmtg(m,labels) 
% -------------------------------------------------------------------------
% Anne-Sophie Dubarry
% 2013/10/01 
%
% DESCRIPTION
%   Function aims to subtract pairs of monopolar SEEG (bipolar montage)
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
% Copyright (C) 2016-2021 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% History : 
% 2013-10-01 :  ASD creation
%
% -------------------------------------------------------------------------
% TODO : display warning message when missing contacts makes important gaps 

% Chaeck number of labels = number of sens in m
if length(labels)~=size(m,1)
    error('Number of momnopolar channels must equal number of labels');
end

% Sort SEEG lables in order (prefix + indices) 
[~,res_index] = mia_sort_nat(labels, 'ascend');
labelo = labels(res_index) ;

a = regexprep(labelo,'\d+','0');    
% idx = ismember(
ua=unique(a) ; 
% 
b = [] ;
idx1 = [] ; idx2 = [] ; 
bilabels={};
% loop over seeg electrodes
for ii=1:length(ua)

  c = find(ismember(a,ua(ii)));
%   bilabels(length(bilabels)+1:length(bilabels)+length(c)-1) = strcat(labels(c(2:end)),'_',labels(c(1:end-1)));
  bilabels(length(bilabels)+1:length(bilabels)+length(c)-1) = strcat(labelo(c(2:end)),'_',labelo(c(1:end-1)));
  idx1 = cat(2,idx1,c(2:end));
  idx2 = cat(2,idx2,c(1:end-1));
  
  %     b = [b; m(c(2:end),:,:) - m(c(1:end-1),:,:)];
  b = [b; m(res_index(c(2:end)),:,:) - m(res_index(c(1:end-1)),:,:)];
  
end
