function scFrames(ax, data)
% scFrames standard sigTOOL line plot function. Single/Multiple Frames
%
% Example:
% scFrames(data)
% scFrames(ax, data)
%   plot data to the current axes, or those specified in ax.
%
% Data should be a data element from a sigTOOLResultData object
% -------------------------------------------------------------------------
% Author: Malcolm Lidierth 01/07
% Copyright ? The Author & King's College London 2007-
% -------------------------------------------------------------------------

if nargin==1
    ax=gca;
    data=ax;
end
view(ax,2);
h1=[];
h2=[];

if size(data.rdata, 2)~=size(data.tdata,2)
%    error('Unmatched vector lengths');
end

for k=1:size(data.rdata, 1);
    
    if isfield(data, 'errdata') && ~isempty(data.errdata.r)
        if isreal(data.errdata.r)
            h1=line('Parent', ax,...
                'XData',data.tdata,...
                'YData', data.rdata(k, :)+data.errdata.r(k,:),...
                'Color', [0.5 0.5 0.5],...
                'Tag', 'sigTOOL:ErrorData',...
                'Visible', 'off',...
                'UserData', k);
            h2=line('Parent', ax,...
                'XData',data.tdata,...
                'YData', data.rdata(k, :)-data.errdata.r(k,:),...
                'Color', [0.5 0.5 0.5],...
                'Tag', 'sigTOOL:ErrorData',...
                'Visible', 'off',...
                'UserData', k);
        else
            h1=line('Parent', ax,...
                'XData',data.tdata,...
                'YData', data.errdata.r.upper(k,:),...
                'Color', [0.5 0.5 0.5],...
                'Tag', 'sigTOOL:ErrorData',...
                'Visible', 'off',...
                'UserData', k);
            h2=line('Parent', ax,...
                'XData',data.tdata,...
                'YData', data.errdata.r.lower(k,:),...
                'Color', [0.5 0.5 0.5],...
                'Tag', 'sigTOOL:ErrorData',...
                'Visible', 'off',...
                'UserData', k);
        end
    end

    if isfield(data, 'barFlag') && data.barFlag==true
        [x,y]=stairs(data.tdata, data.rdata(k, :));
        h3=line('Parent', ax,...
        'XData', x,...
        'YData', y,...
        'Color', [0.1 0.1 0.5],...
        'Tag', 'sigTOOL:ResultData',...
        'UserData', k);
    else
    h3=line('Parent', ax,...
        'XData',data.tdata,...
        'YData', data.rdata(k, :),...
        'Color', [0.1 0.1 0.5],...
        'Tag', 'sigTOOL:ResultData',...
        'UserData', k);
    end
    
    grp=hggroup('Parent',ax);
    set([h1 h2 h3], 'Parent', grp);
    set(grp,'Tag','sigTOOL:ResultGroup');
end

if size(data.rdata, 1)==1
    h=findobj(ax, 'Tag', 'sigTOOL:ErrorData');
    if ~isempty(h)
        set(h, 'Visible', 'on');
    end
end

setappdata(ancestor(ax,'figure'),'sigTOOLViewStyle', '2D')
return
end