function out=horzcat(varargin)
% HORZCAT method overloaded for adcarray objects
%
% Example:
% C = vertcat(A1, A2, ...)
% returns a double precision array C. 
% Called by matlab when A1, A2 etc includes an adcarray object
%
% See also HORZCAT, ADCARRAY/VERTCAT
%
% Author: Malcolm Lidierth
% Copyright ? The Author & King's College London 2006

index(1).type='()';
index(1).subs={};

out=builtin('horzcat',subsref(varargin{1},index),subsref(varargin{2},index));
if nargin>2
    for i=3:nargin
        out=builtin('horzcat',out,subsref(varargin{i},index));
    end;
end;
