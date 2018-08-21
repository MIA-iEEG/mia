function [h,offset] = stackplot_ASD(t,x,color,labels, offset)
%STACKPLOT Plot in a stacked fashion
% function [h,offset] = stackplot(t,x,offset)
% t is time (optional if offset NOT given)
% x is the data matrix, one plot per column
% offset is the desired offset, optional, useful primarily to align two
% sets of data, i.e.
%  [h,offset] = stackplot(t,x);
%  hold on, stackplot(t,x2,offset); hold off
%
% output is optional, returning handles and offset
% Uses: stackplot(x), stackplot(t,x), stackplot(t,x,offset)
%  NOT: stackplot(x,offset) (cbb)

% -------------- Script History --------------
% 1990s original versions by JCM
% Sep 2006 JCM redone again
% 4-Oct-2006 JCM Commenting, designing different separations
% Jan 2013 ASD adding color + labels 
% Jan 2016 : ASD Flip display order (-offset and fliplr(labels)
% --------------------------------------------



if nargin == 2,
    color = x ;
	x = t;
	t = (1:size(x,1))';
    
end


if nargin < 5,
    offset = []; % not given
end

[N,PLOTS] = size(x); % length of data, number of plots

if isempty(offset), % need to calculate an offset
    % two flavors of calculating sig
    switch 'each'
        case 'overall'
            % calculate an overall standard dev applied to all plots
            SEP_SIGS = 4; % number of standard deviations to separate plots
            sig = std(x(:)); % standard dev of all the data
            sig = SEP_SIGS*(0:(PLOTS-1))*sig; % same offset
        case 'each'
            % each plot gets a unique std. dev.
            SEP_SIGS = 2; % number of standard deviations to separate plots
            sig = std(x,0,1); % standard dev of each column
            if PLOTS > 1
                sig = sum(hankel(sig(1:2),sig(2:end)))*SEP_SIGS;
                sig = [0 cumsum(sig)];
            end

    end
    offset = repmat(sig,N,1);
end % else use what user gave

if isreal(x),
    hf = plot(t,x - offset,'LineWidth',1.05);
%     hf = plot(t,x + offset,color);
    % Add grid lines
    HOLDON = ishold; % is hold on or off already?
    if ~HOLDON
        hold on
    end
    % only need beginning and ending points for efficiency
    plot(t([1 end]),-offset([1 end],:),':','Color','k')
    
    if ~HOLDON
        % hold was not on already
        hold off
    end
    
else
    % complex
    hf = plot(x + sqrt(-1)*offset); % vertically stacked
end


axis tight

if ~isempty(labels)
    set(gca,'YTick',fliplr(-1*offset(1,:)),'YTickLabel',fliplr(labels)) ;
end

if nargout, % user wants handles
    h = hf;
end
