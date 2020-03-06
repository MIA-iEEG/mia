function [seg,r_edge] = mia_get_significant_segments(h) 
% ***********************************************************************
%  [] = get_significant_segments(h) 
% Returns the lengths of consecutive significant samples in h
% 
% ***********************************************************************
%
%
% ***********************************************************************
% Inputs :
% h : nsens x nsamples 
%
% Output :
% seg : durations of significant segments (in samples)
% r_edge : Position of the segments in a vector shape
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
% Copyright (C) 2016-2020 CNRS - Universite Aix-Marseille
%
% ========================================================================
% This software was developed by
%       Anne-Sophie Dubarry (CNRS Universite Aix-Marseille)
% 2014/10/29 : ASD : optimization, w vectors
% 

seg = [] ;
r_edge = [];

% Add 0 at the edges of h (if any significant window not closed) 
h(:,1) = 0 ; h(:,end)=0;

% Create on vector 
vc = reshape(h',1,[]) ;

% Finds rising and falling edges 
r_edge = find(diff(vc')'>0) ;
f_edge = find(diff(vc')'<0) ;

% Segments durations
seg = f_edge - r_edge ;
