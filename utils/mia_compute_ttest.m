function [tarray, parray] = mia_compute_ttest(dc) 
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
% Modified by JCM October 2015 to not require STATS TTEST, but rather run
% the local function attached.
% -------------------------------------------------------------------------

% Initialize stat arrays (tvals and pvals)
tarray = zeros(size(dc,1),size(dc,2));
parray = zeros(size(dc,1),size(dc,2));


% Loop over contact
for contactidx=1:size(dc,1)

    %[~,P,~,STATS] = ttest(squeeze(dc(contactidx,:,:))');
    [P,T] = ttest(squeeze(dc(contactidx,:,:))'); % JCM
    %tarray(contactidx,:) = STATS.tstat ; 
    tarray(contactidx,:) = T ; % JCM
    parray(contactidx,:) = P; 
end


function [p,t] = ttest(x) % JCM

% Make sure there is no NaN otherwise tvalues = 0 and pvalues =1 
if  ~sum(isnan(x(:)))
    t = mean(x) ./ sqrt(std(x).^2 / size(x,1));
    n = size(x,1)-1;
    p = betainc(n./(n+t.^2),0.5*n,0.5);
else
    t = zeros(1,size(x,2));
    p = ones(1,size(x,2));
end
