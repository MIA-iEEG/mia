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
function [] = display_raster(varargin)

% ASD 2015/5/119

if nargin == 5
    t=varargin{1};
    mat=varargin{2};
    rt=varargin{3};
    colorpt=varargin{4};
    ax=varargin{5};
elseif nargin==3;
    t=varargin{1};
    mat=varargin{2};
    ax=varargin{3};
end

h = filter_timewin_signif(abs(mat')>0, 10) ;

m = mat'.*h ;

% Filters the tvals
for jj=1:size(m,2) ; vec(:,jj) = filtfilt(ones(1,2)/2,1,m(:,jj)) ; end ;

axes(ax);
imagesc(t,1:size(mat,2),vec) ; colorbar ; caxis([-10 10]) ;
hold on ;

if nargin==5
    for pp=1:length(rt)
        plot(rt(pp),pp,'.','Color',colorpt(pp,:),'MarkerSize',20);
    end
end
set(gca,'XGrid','on') ; %title(gca,labroi, 'FontSize', 15);
set(gcf,'color','white')    ;  xlabel(gca,'Time (s)') ; set(gca,'Fontsize',10);
    

end
