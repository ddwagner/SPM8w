function [hours,minutes,seconds] = spm8w_timecalc(TimeInSeconds)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: December, 2009 - DDW
% ==============================================================================
% FORMAT [hours,minutes,seconds] = spm8w_timecalc(TimeInSeconds);
% 
% Simple function to return time in hours, minutes and seconds given input
% in seconds. Useful for determining how long stages in processing lasted
% in human readable terms.
% ==============================================================================
% CHANGE LOG:
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Convert seconds to HH:MM:SS
seconds = TimeInSeconds;
hours = fix(seconds/3600);
seconds = seconds - 3600*hours; % get number of hours
minutes = fix(seconds/60);      % get number of minutes
seconds = seconds - 60*minutes; % remove the minutes
seconds = round(seconds);