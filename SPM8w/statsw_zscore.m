function [z] = statsw_zscore(x)
% ==============================================================================
% STATSW
% Quick and dirty replacements for functions in matlab stats toolbox required
% by various functions in SPM8w. 
% 
% Heatherton & Kelley Labs
% Last update: March, 2013 - DDW
% Created: March, 2013 - DDW
% ==============================================================================
% statsw_zscore(data)
%
% Returns the z-score of the input data (i.e., data-mean(data) ./std(data)).
% The resultant z-score vector will be of same length as data.
% ==============================================================================
% CHANGE LOG:
% -First version - March, 2013
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Input checks
switch (nargin)
  case nargin>1
    error('Please provide only a single vector of data');
end

%---Check format (we want row vectors)
% if size(x,2) < size(x,1)
%     x = x';
% end

%---Zscore is data-mean(data) ./std(data)
z = (x-mean(x)) ./ std(x);