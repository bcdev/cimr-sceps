function y = nanmedian(varargin)
%NANMEDIAN Median value, ignoring NaNs.
%   M = NANMEDIAN(X) returns the sample median of X, treating NaNs as missing
%   values.  For vector input, M is the median value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the median value of
%   non-NaN elements in each column.  For N-D arrays, NANMEDIAN operates
%   along the first non-singleton dimension.
%
%   NANMEDIAN(X,DIM) takes the median along dimension DIM of X.
%
%   See also MEDIAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2016 The MathWorks, Inc.


narginchk(1,2);
y = median(varargin{:},'omitnan');
