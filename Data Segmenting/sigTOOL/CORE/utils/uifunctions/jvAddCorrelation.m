function h=jvAddCorrelation(h)
% jvAddCorrelation addpanel function
% 
% Example:
% h=jvAddCorrelation(h)
%
% -------------------------------------------------------------------------
% Author: Malcolm Lidierth 11/07
% Copyright ? The Author & King's College London 2007-
% -------------------------------------------------------------------------


Height=0.09;
Top=0.75;

h=jvAddPanel(h, 'Title', 'Options',...
    'dimension', 0.5);

h=jvElement(h{end},'Component', 'javax.swing.JComboBox',...
    'Position',[0.1 Top 0.8 Height],...
    'DisplayList', {'0.1' '0.2' '0.5' '1.0' '2.0'},...
    'Label', 'Maximum Lag (s)',...
    'ToolTipText', 'Data sections length (ms)');

h=jvElement(h{end},'Component', 'javax.swing.JComboBox',...
    'Position',[0.1 Top-(2*Height) 0.8 Height],...
    'DisplayList', {  'Coeff' 'Unbiased' 'Biased' 'None'},...
    'Label', 'Scaling',...
    'ToolTipText', 'Result scaling mode');

h=jvElement(h{end},'Component', 'javax.swing.JCheckBox',...
    'Position',[0.1 Top-(3.5*Height) 0.8 Height],...
    'Label', 'Remove DC',...
    'ToolTipText', 'Pre-process: Remove DC & linear trend');

h=jvElement(h{end},'Component', 'javax.swing.JComboBox',...
    'Position',[0.1 Top-(6*Height) 0.8 Height],...
    'DisplayList', {'22' '21' '20' '19' '18'},...
    'Label', 'Max BlockSize(2^x)',...
    'ToolTipText', 'Data block size');

return
end
