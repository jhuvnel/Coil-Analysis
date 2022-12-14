function chantitle=SONGetChanTitle(fh, chan)
% SONGETCHANTITLE Returns the channel title
%
% TITLE=SONGETCHANTITLE(FH, CHAN)
% 
% INPUTS FH is the file handle
%         CHAN is the channel number (0 to SONMaxChans()-1
% OUTPUT the channel title as a string
%
% Author:Malcolm Lidierth
% Matlab SON library:
% Copyright ? The Author & King's College London 2005-2006
             
chantitle='0123456789';
chantitle=calllib('son32','SONGetChanTitle', fh, chan, chantitle);
return;