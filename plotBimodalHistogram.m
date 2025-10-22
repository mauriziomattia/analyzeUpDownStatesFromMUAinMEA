function ModeParams = plotBimodalHistogram(Values, ValRange, BinNum)
%
%  ModeParams = plotBimodalHistogram(Values[, ValRange[, BinNum]])
%
%
%   Copyright 2009 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Sep. 22, 2009
%

GRAY_CLR = repmat(0.9, 1, 3);
BIN_NUM = 50;
PEAK_THRESHOLD = 3/5; % The level from the peak which is at about 1 st.dev from it.
TAIL_TO_NEGLECT = 0.005;
MIN_2ND_PEAK_AREA = 10.0;

%
% Parameter settings...
SampleSize = length(Values);

if exist('BinNum')
   if BinNum>0
      BIN_NUM = BinNum;
   end
end

if exist('ValRange') ~= 1
   ValRange = [min(Values) max(Values)];
else
   if isempty(ValRange)
      ValRange = [min(Values) max(Values)];
   end
end
X = linspace(ValRange(1), ValRange(2), BIN_NUM);
Values = Values(find(Values >= ValRange(1) & Values <= ValRange(2)));

%
% Computes histogram...
N = hist(Values, X) / SampleSize * 100;

%
% Smooth the histogram in order to find maxima...
knt = augknt(linspace(ValRange(1), ValRange(2), round(BIN_NUM/4)), 4);
spN = spap2(knt, 4, X, N);
spN = spap2(newknt(spN), 4, X, N);
spN = spap2(newknt(spN), 4, X, N); % A second refinement of the fit...

% removes tails of the distribution for low and high values...
Y = cumsum(N);
Y = Y / Y(end);
a = X(find(Y > TAIL_TO_NEGLECT, 1, 'first'));
b = X(find(Y <= 1.0 - TAIL_TO_NEGLECT, 1, 'last'));

% finds extrema (minima and maxima) in the bulk of the distribution...
XMaxMin = mean(fnzeros(fnder(spN), [a b]));
YMaxMin = fnval(spN, XMaxMin);

% Criterium 1 for peak selection (original):
[YMax, YMaxNdx] = max(YMaxMin);

% Criterium 2 for peak selection (18th Sep. 2009):
% ConcMaxMin = fnval(fnder(fnder(spN)),XMaxMin);
% [ConcMax, YMaxNdx] = min(ConcMaxMin);

%
% Gaussian fit on the first peak...
XFirstMax = XMaxMin(YMaxNdx);
YFirstMax = YMaxMin(YMaxNdx);
if XFirstMax-a < b-XFirstMax  % Is on left or on the right?
   FirstLeft = 1;
   ndx = find(X >= XFirstMax, 1 , 'first');
   [val, pos] = max(N(ndx-1:ndx+1));
   ndx = ndx + pos - 2;
   ndxL = find(N(1:ndx)>val*PEAK_THRESHOLD, 1, 'first');
%    ndxL = find(N(1:ndx)>YFirstMax*PEAK_THRESHOLD, 1, 'first');
   
   beta0(1) = X(ndx);
   beta0(2) = (ndx - ndxL)*(X(2) - X(1));
%    if beta0(2) == 0
%       beta0(2) = X(2) - X(1);
%    end
   beta0(3) = val;
   beta0(3) = YFirstMax;
   % Original
   ndx = ndx + (ndx - ndxL - 1);
   % 18th Sep. 2009
%    ndx = ndx + (ndx - ndxL - 1) + 1;
   beta1 = nlinfit(X(1:ndx), N(1:ndx), @mygauss, beta0);

   ModeParams.Mu1 = beta1(1);
   ModeParams.Sigma1 = beta1(2);
   ModeParams.Ampl1 = beta1(3);
else
   FirstLeft = 0;
   ndx = find(X <= XFirstMax, 1 , 'last');
   [val, pos] = max(N(ndx-1:ndx+1));
   ndx = ndx + pos - 2;
   ndxL = find(N(ndx:end)<val*PEAK_THRESHOLD, 1, 'first');
%    ndxL = find(N(ndx:end)<YFirstMax*PEAK_THRESHOLD, 1, 'first');

   beta0(1) = X(ndx);
   beta0(2) = ndxL*(X(2) - X(1));
%    if beta0(2) == 0
%       beta0(2) = X(2) - X(1);
%    end
   beta0(3) = val;
   beta0(3) = YFirstMax;
   % Original
  ndx = ndx - ndxL;
   % 18th Sep. 2009
%    ndx = ndx - ndxL - 1;
   beta2 = nlinfit(X(ndx:end), N(ndx:end), @mygauss, beta0);

   ModeParams.Mu2 = beta2(1);
   ModeParams.Sigma2 = beta2(2);
   ModeParams.Ampl2 = beta2(3);
end

%   
% Gaussian fit on the second peak...
if FirstLeft == 0  % Is on left or on the right?
   N2 = N - mygauss(beta2, X);
   SecondPeakArea = sum(N2(find(N2>0)));
   disp(['2nd peak area = ' num2str(SecondPeakArea) '%']);
   if SecondPeakArea >= MIN_2ND_PEAK_AREA
      % Criterium 0 for peak selection:
      %    [val, ndx] = max(N2);

      % Criterium 1 for peak selection (original):
      ndx = find(X>a);
      maxndx = find(diff(N2(ndx))<0, 1, 'first');
      ndx = ndx(maxndx);
      val = N2(ndx);

      % Criterium 2 for peak selection (18th Sep. 2009):
%       ndx = find(X>=N2*X'/sum(N2), 1, 'first');
%       val = N2(ndx);

      ndxL = find(N2(1:ndx)>val*PEAK_THRESHOLD, 1, 'first');

      % Criterium 3 for peak selection (18th Sep. 2009):
%       spN2 = fnval(spN, X) - mygauss(beta2, X);
%       [val, ndx] = max(spN2);
%       
%       ndxL = find(spN2(1:ndx)>val*PEAK_THRESHOLD, 1, 'first');

      beta0(1) = X(ndx);
      beta0(2) = (ndx - ndxL)*(X(2) - X(1));
      beta0(3) = val;
      ndx = ndx + (ndx - ndxL - 1);
%       ndx = ndx + (ndx - ndxL - 1) + 1;
      ndx = min([ndx length(X)]);
      beta1 = nlinfit(X(1:ndx), N2(1:ndx), @mygauss, beta0);

      ModeParams.Mu1 = beta1(1);
      ModeParams.Sigma1 = beta1(2);
      ModeParams.Ampl1 = beta1(3);
   end
else
   N2 = N - mygauss(beta1, X);
   SecondPeakArea = sum(N2(find(N2>0)));
   disp(['2nd peak area = ' num2str(SecondPeakArea) '%']);
   if SecondPeakArea >= MIN_2ND_PEAK_AREA
      % Criterium 0 for peak selection:
      %    [val, ndx] = max(N2);

      % Criterium 1 for peak selection (original):
      ndx = fliplr(find(X<b));
      maxndx = find(diff(N2(ndx))<0, 1, 'first');
      ndx = ndx(maxndx);
      val = N2(ndx);

      % Criterium 2 for peak selection (18th Sep. 2009):
%       ndx = find(X>=N2*X'/sum(N2), 1, 'first');
%       val = N2(ndx);

      ndxL = find(N2(ndx:end)<val*PEAK_THRESHOLD, 1, 'first');

      % Criterium 3 for peak selection (18th Sep. 2009):
%       spN2 = fnval(spN, X) - mygauss(beta1, X);
%       [val, ndx] = max(spN2);
%       
%       ndxL = find(spN2(ndx:end)<val*PEAK_THRESHOLD, 1, 'first');

      beta0(1) = X(ndx);
      beta0(2) = ndxL*(X(2) - X(1));
      beta0(3) = val;
      ndx = ndx - ndxL - 1;
      beta2 = nlinfit(X(ndx:end), N2(ndx:end), @mygauss, beta0);

      ModeParams.Mu2 = beta2(1);
      ModeParams.Sigma2 = beta2(2);
      ModeParams.Ampl2 = beta2(3);
   end
end

%
% Plots the histogram...
figure; 
hold on;

[XX, YY] = stairs(X - (X(2) - X(1))/2, N);
XX = [XX(1) XX' XX(end)];
YY = [0 YY' 0];
patch(XX, YY, GRAY_CLR, 'FaceColor', GRAY_CLR, 'EdgeColor', GRAY_CLR);
plot(XX, YY, 'k', 'LineWidth', 1);

FitHist = zeros(size(N));
if exist('beta1')
   plot(X, mygauss(beta1, X), 'b', 'LineWidth', 1);
   FitHist = mygauss(beta1, X);
end
if exist('beta2')
   plot(X, mygauss(beta2, X), 'r', 'LineWidth', 1);
   FitHist = FitHist + mygauss(beta2, X);
end
plot(X, FitHist, 'g', 'LineWidth', 1);

YLim = get(gca, 'YLim');

if exist('beta1')
   plot(ModeParams.Mu1 + [0 0], [0 YLim(2)], 'b--');
   text(ModeParams.Mu1, ModeParams.Ampl1+0.5, ...
        [' ' num2str(ModeParams.Mu1, 3) ...
         ' (s.d.: ' num2str(ModeParams.Sigma1, 3) ')']);
end
if exist('beta2')
   plot(ModeParams.Mu2 + [0 0], [0 YLim(2)], 'r--');
   text(ModeParams.Mu2, ModeParams.Ampl2+0.5, ...
        [' ' num2str(ModeParams.Mu2, 3) ...
         ' (s.d.: ' num2str(ModeParams.Sigma2, 3) ')']);
end

set(gca, 'Layer', 'top', 'Box', 'off', 'TickDir', 'out');
set(gca, 'XLim', ValRange);

xlabel('Value (a.u)');
ylabel('Samples (%)');


%
% Work out the error of the fit...
%
HistError = N/100.;
HistError = sqrt((HistError .* (1. - HistError)) / length(Values));
ndx = find(abs(N-FitHist)/100. > HistError);
ModeParams.FitError = sum(abs(N(ndx)-FitHist(ndx))/100. - HistError(ndx));

