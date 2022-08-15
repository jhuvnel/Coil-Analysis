function varargout=menu_CopyObjects(varargin)
% menu_CopyObjects helper for copying graphics objects
% 
% menu_CopyObjects(hObject, EventData)
%     standard menu callback
%     
% -------------------------------------------------------------------------
% Author: Malcolm Lidierth 11/06
% Copyright � The Author & King's College London 2006-
% -------------------------------------------------------------------------
if nargin==1 && varargin{1}==0
    varargout{1}=true;
    varargout{2}='Copy MATLAB Objects';
    varargout{3}=0;
    return
end

CopyObjects(gca);
