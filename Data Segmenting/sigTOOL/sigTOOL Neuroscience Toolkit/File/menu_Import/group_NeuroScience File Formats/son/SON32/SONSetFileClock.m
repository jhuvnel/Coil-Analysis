function SONSetFileClock(fh, usPerTime, timePerADC)
% SONSETFILECLOCK sets the basic time units and the clocks per ADC conversion
% 
% SONSETFILECLOCK(FH, USPERTIME, TIMEPERADC)
% INPUTS: FH the SON file handle
%         USPERTIME the number of microseconds per clock tick
%         TIMEPERADC the number of clock ticks per ADC conversion
%         
% See CED documentation for further details
% 
% Author:Malcolm Lidierth
% Matlab SON library:
% Copyright ? The Author & King's College London 2005-2006

calllib('son32','SONSetFileClock', fh, usPerTime, timePerADC);
