function axes_equal_loglog ( hAxes,xLimits,yLimits )
%This Function replaces axes equal for loglog plot as the former doesn't
%present acceptaveble results when used with the latter.
%If you don't want to set the axis limits, input only hAxis.
%Example : axes_equal_loglog (gca) where gca the current axis handle

% In case I don't want to limit the axis
if (nargin < 2) || isempty(xLimits)
    xLimits = get(hAxes,'XLim');
    msg = 'xLimits may be empty';
    warning(msg);
end
if (nargin < 3) || isempty(yLimits)
    yLimits = get(hAxes,'YLim');
    msg = 'yLimits may be empty ';
    warning(msg);
end

% We need to consider both number and sizes of decades in the axis range
% Therefore, we calculate the average size for the decade of each axis and
% then using the ration of that average for both axis.

% i.e. diff(xLimits) total size of the x axis and diff(log10(xLimits) the
% number of decades

averagey = diff(yLimits)/diff(log10(yLimits));
averagex = diff(xLimits)/diff(log10(xLimits));

set(hAxes,'XLim',xLimits,'YLim',yLimits,'DataAspectRatio',[1 averagey/averagex 1]);


end

