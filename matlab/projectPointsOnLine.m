function [xProj, yProj, varargout] = projectPointsOnLine(x, y, xLine, yLine)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                              %
% function [xProj, yProj, varargout] = projectPointsOnLine(x, y, xLine, yLine) %
%                                                                              %
% This function projects the given points on a fitted line.                    %
% [xProj, yProj] = projectPointsOnLine(x, y, xLine, yLine).                    %
% fits a straight line with the points defined in (xLine, yLine)               %
% and finds the projection of the points (x, y) over the fitted line.          %
% [xProj, yProj, xFitLine, yFitLine] = projectPointsOnLine(x, y, xLine, yLine).%
% returns also the line fit from the original set of points                    %
%                                                                              %
% Author: Ananda Pascual (IMEDEA), SOCIB team (www.socib.es)                   %
%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Adjust a straight line from the initial points
% (user could provide more than two and not in perfect straight line)
% Write the line equation y = mx + b in matrix form y = [x, 1] * [m; b]
% and finds the line coefficients (slope and intercept)
xMat = [xLine(:), ones(numel(xLine), 1)];
lineCoefs = xMat \ yLine(:);
yLine = xMat * lineCoefs;

% If user asked for it, return also the fitted line
if nargout > 1,
    varargout{1} = xLine;
    varargout{2} = yLine;
end;

% Get first and last points, so the line will be defined by two endpoints
xLine = [xLine(1); xLine(end)];
yLine = [yLine(1); yLine(end)];

% Find a vector that defines the line direction
line = [xLine(:), yLine(:)];
lineVector = diff(line, 1, 1);
unitLineVector = lineVector ./ sqrt(sum(lineVector.^2));

% Find a vector that defines how the point departs from the line
outerVector = [x - xLine(1), y - yLine(1)];

% Compute the dot product between the line vector and the outer vector 
% to find the projection of the outer vector over the line
repeatLineVector = ones(length(outerVector),1) * unitLineVector;
dist = sum(repeatLineVector .* outerVector, 2);

xProj = xLine(1) + dist * unitLineVector(1);
yProj = yLine(1) + dist * unitLineVector(2);

return;
