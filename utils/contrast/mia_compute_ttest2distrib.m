function [tarray, parray] = mia_compute_ttest2distrib(dc1,dc2) 
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
% -------------------------------------------------------------------------
% DESCRIPTION
%   Perform t-test against 0 
%
% Inputs :
%           d : Matrix containing the data
%           (sensors x time samples x trials)
% 
% Output:   t-values and p-value
% 
% -------------------------------------------------------------------------

% Initialize stat arrays (tvals and pvals)
tarray = zeros(size(dc1,1),size(dc1,2));
parray = zeros(size(dc1,1),size(dc1,2));

% Loop over contact
for contactidx=1:size(dc1,1)

    [tstat df] = mia_ttest2_cell({ squeeze(dc1(contactidx,:,:)) squeeze(dc2(contactidx,:,:)) });
 
%     [~,P,~,STATS] = ttest2(squeeze(dc1(contactidx,:,:))',squeeze(dc2(contactidx,:,:))');
    tarray(contactidx,:) = tstat ; 
    parray(contactidx,:) = 1-mia_tcdf(tstat,df-1) ;% Probability of larger t-statistic    
%     end
end
