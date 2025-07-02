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
% -----------------------------------------------------------------------
function [single_rois] = mia_n_to_one_roi(rois,Conditions) 

% Init (fileds will be replaced) 
single_rois = rois(1,:);

% Preallocate or initialize a list of matrices
allMatrices = cellfun(@(s) s.signmoy, rois, 'UniformOutput', false);  % extract all matrices
avgMatrices = cellfun(@(a,b) (a + b) / 2, allMatrices(1,:), allMatrices(2,:), 'UniformOutput', false);
single_rois = cellfun(@(s, d) setfield(s, 'signmoy', d), single_rois , avgMatrices, 'UniformOutput', false);

% Preallocate or initialize a list of matrices
allMatrices = cellfun(@(s) s.Fmask, rois, 'UniformOutput', false);  % extract all matrices
avgMatrices = cellfun(@(a,b) (a + b) / 2, allMatrices(1,:), allMatrices(2,:), 'UniformOutput', false);
single_rois= cellfun(@(s, d) setfield(s, 'Fmask', d), single_rois , avgMatrices, 'UniformOutput', false);

% Preallocate or initialize a list of matrices
allMatrices = cellfun(@(s) s.F, rois, 'UniformOutput', false);  % extract all matrices
avgMatrices = cellfun(@(a,b) (a + b) / 2, allMatrices(1,:), allMatrices(2,:), 'UniformOutput', false);
single_rois= cellfun(@(s, d) setfield(s, 'F', d), single_rois , avgMatrices, 'UniformOutput', false);

% Extract the signmoy, Fmask, and F fields from each ROI into separate cell arrays
allSignmoyCell = cellfun(@(s) s.signmoy, rois, 'UniformOutput', false);
allFmaskCell = cellfun(@(s) s.Fmask, rois, 'UniformOutput', false);
allFCell = cellfun(@(s) s.F, rois, 'UniformOutput', false);

for iRois=1:size(rois,2)

    single_rois{iRois}.name = sprintf('%s - ',single_rois{iRois}.name,Conditions) ;
      single_rois{iRois}.signmoyAll = allSignmoyCell(:,iRois);  % Assign signmoyAll field
    single_rois{iRois}.FmaskAll = allFmaskCell(:,iRois);  % Assign FmaskAll field
    single_rois{iRois}.FAll = allFCell(:,iRois);  % Assign FAll field
end

