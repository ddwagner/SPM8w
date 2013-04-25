function spm8w_timeplot(varargin)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: April, 2012 - DDW
% ==============================================================================
% spm8w_timeplot(tseries, name, TR, tylim, fylim, printfig)
% 
% Simple function to make figures of timeseries and their power spectra.
% Takes a timeseries vector, its name (i.e.'Amygdala' or 'Linear Trend',
% etc.) and the TR of the experiment.
%
% Optionally takes YLimits for the Time and Frequency domain plots. 
% Good settings are [-4,4] for Time domain and [0,0.04] for the Frequency 
% domain. 
%
% Example: 
% spm8w_timeplot(LIN,'Linear', 2.2)
% spm8w_timeplot(POLY2,'Filtering with Poly2', 2.2,[-5,5],[0,0.1]);
% ==============================================================================
% CHANGE LOG:
%
% =======1=========2=========3=========4=========5=========6=========7=========8

%---Input checks
switch (nargin)
  case nargin < 3
    error('Too few paramters. You must specify the timeseries, name and TR');
  case 3
    tseries = varargin{1};
    name    = varargin{2};
    TR      = varargin{3};
    tylim   = [];
    fylim   = [];
    prntfig = 0;
  case 4
    tseries = varargin{1};
    name    = varargin{2};
    TR      = varargin{3};
    tylim   = varargin{4};
    fylim   = [];
    prntfig = 0;
  case 5    
    tseries = varargin{1};
    name    = varargin{2};
    TR      = varargin{3};
    tylim   = varargin{4};
    fylim   = varargin{5};
    prntfig = 0;
  case 6    
    tseries = varargin{1};
    name    = varargin{2};
    TR      = varargin{3};
    tylim   = varargin{4};
    fylim   = varargin{5};
    prntfig = varargin{6};
  otherwise 
    error('Too many paramters.'); 
end

%---Plot Timeseries
figure1 = figure('Units','pixels','Position',[100, 100, 800, 600]);
subplot(2,1,1);
plot(tseries,'LineWidth',1,'Color','r');
xlabel(['Scans TR:',num2str(TR)])
ylabel('Timeseries')
title('Time Domain','Interpreter','Tex');
grid on
axis tight
if ~isempty(tylim)
    ylim(tylim);
end
legend(name)

%---Plot Timeseries
subplot(2,1,2)
HPF   = 128;
%need to meancenter for fft
if abs(mean(tseries)) > 1
    tseries = spm_detrend(tseries);
end
gX    = abs(fft(tseries)).^2;
gX    = gX*diag(1./sum(gX));
q     = size(gX,1);
Hz    = [0:(q - 1)]/(q*TR);
q     = 2:fix(q/2);
plot(Hz(q),gX(q,:),'r')
if isempty(fylim)
    patch([0 1 1 0]/HPF,[0 0 1 1]*max(max(gX)),[1 0.4 0.4]*.3,'facealpha',.5);
else
    patch([0 1 1 0]/HPF,[0 0 1 1]*fylim(2),[1 0.4 0.4]*.3,'facealpha',.5);   
end
xlabel('Frequency (Hz)')
ylabel('relative spectral density')
title('Frequency domain','Interpreter','Tex');
legend(name)
grid on
axis tight
if ~isempty(fylim)
    ylim(fylim);
end

%---Save Figure
if prntfig
    set(figure1, 'PaperPositionMode', 'auto');
    print(figure1, '-dpng','-r120',[name,'.png']);
end
