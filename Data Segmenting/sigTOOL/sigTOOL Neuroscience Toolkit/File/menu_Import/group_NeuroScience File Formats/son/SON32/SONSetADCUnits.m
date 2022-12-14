function SONSetADCUnits(fh, chan, units)
% SONSETADCUNITS sets the units string for a channel
% 
% SONADCUNITS(FH, CHAN, UNITS)
% where FH is the SON file handle, CHAN is the channel number and UNITS is
% the string  to be written to the file.
% 
% No return value
%
% Author:Malcolm Lidierth
% Matlab SON library:
% Copyright ? The Author & King's College London 2005-2006

global SON_UNITSZ;

units=[units(1:min(length(units),SON_UNITSZ)) ''];
size(units)
calllib('son32','SONSetADCUnits', fh, chan, units);