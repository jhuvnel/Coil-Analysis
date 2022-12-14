function SONSetFileComment(fh, which, comment)
% SONSETFILECOMMENT sets one of the five file comment fields
% 
% SONSETFILECOMMENT(FH, WHICH, COMMENT)
%         FH is the SON file handle
%         WHICH selects the comment (0 - 4)
%         COMMENT is a string
% 
% No return value
%
% Author:Malcolm Lidierth
% Matlab SON library:
% Copyright ? The Author & King's College London 2005-2006

global SON_COMMENTSZ;

comment=comment(1:min(SON_COMMENTSZ,length(comment)));
calllib('son32','SONSetFileComment', fh, which, comment);
