% [data, npoints, s]=ABFGetDACChannel(dllname, filename, channel, episode)
% ABFGetDACChannel reads a data channel from an Axon Instruments ABF file
% 
% Examples:
% [data, npoints, s]=ABFGetDACChannel(dllname, filename, channel);
% [data, npoints, s]=ABFGetDACChannel(dllname, filename, channel, episode);
% 
% Inputs: 
% dllname: the full path and name of the Axon Instruments ABFFIO.DLL
% filename: the path and file name of the ABF file to read
% channel: the target channel, range 0-15
% episode: if specified the data episode to read
%     
% outputs: the data in a NxM matrix, each column representing one episode
% npoints: a vector, one value for the number of samples in each episode
%           (i.e. each column of data)
% s: a structure with the channel comment, title and units
% 
% data is single floating point and is scales to the units in s.Units
% 
% This routine calls the Axon Instruments DLL and is therefore Windows only
%     
%--------------------------------------------------------------------------
% Author: Malcolm Lidierth 06/07
% Copyright ? The Author & King's College London 2007
%
% Acknowledgements:
% Revisions:
%--------------------------------------------------------------------------