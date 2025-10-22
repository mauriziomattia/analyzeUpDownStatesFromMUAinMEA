function WFs = extractEventCenteredWaveForms(t, Y, Triggers, ...
                                             WFTimeRange, BLTimeRange)
%
%
%   Copyright 2008 Maurizio Mattia @ Ist. Super. Sanità, Rome - Italy
%   Version: 1.0 - Jan. 21, 2008
%


SamplingFreq = 1 / (t(2) - t(1));
ndxRel = round(WFTimeRange(1)*SamplingFreq):round(WFTimeRange(2)*SamplingFreq);
%ndxSurge = round(-DELTA_T/2*SAMPLING_FREQ):round(DELTA_T/2*SAMPLING_FREQ);
if exist('BLTimeRange')
   ndxBaseline = round(BLTimeRange(1)*SamplingFreq):round(BLTimeRange(2)*SamplingFreq);
else
   ndxBaseline = [];
end
tWFs = ndxRel / SamplingFreq;

for n = 1:length(Triggers)
   ndxOffset = round((Triggers(n) - t(1)) * SamplingFreq);
   
   ndx = ndxRel + ndxOffset;
   ndxMask = find(ndx>0 & ndx<=length(t));
   
   WFs(n).time = tWFs(ndxMask);
%   WFs(n).value = csY(ndx(ndxMask)) - mean(csY(ndxSurge + ndxOffset));
%   WFs(n).value = csY(ndx(ndxMask)) - mean(csY(ndxBaseline + ndxOffset));
   if isempty(ndxBaseline)
      WFs(n).value = Y(ndx(ndxMask))';
   else
      WFs(n).value = Y(ndx(ndxMask))' - mean(Y(ndxBaseline + ndxOffset));
   end
   
end
