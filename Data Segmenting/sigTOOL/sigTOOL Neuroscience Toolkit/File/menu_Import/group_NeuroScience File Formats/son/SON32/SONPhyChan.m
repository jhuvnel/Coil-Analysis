function pchan=SONPhyChan(fh, chan)
% SONPHYCHAN returns the physical channel (hardware port) for a channel
% in a file 
% 
% PCHAN=SONPHYCHAN(FH, CHAN)
% where FH is the SON file handle and CHAN is the file channel number.
% PCHAN is negative for channels with no hardware port
% 
% Author:Malcolm Lidierth
% Matlab SON library:
% Copyright ? The Author & King's College London 2005-2006

pchan=calllib('son32','SONPhyChan', fh, chan);