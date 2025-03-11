function mustBeSingleStream(S)

%mustBeSingleStream Validate that STREAMobj contains a single river only
%
% Syntax
%
%     mustBeSingleStream(S)
%
% Description
%
%     The function validates that a STREAMobj contains only a single
%     stream. E.g. trunk(klargestconncomps(S)) will be a single stream.
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 7. March, 2025

if isa(S,'PPS')
    S = S.S;
end

TF = isa(S,'STREAMobj') && info(S,'nrchan')==1;
if ~TF
    error('TopoToolbox:incorrectinput',...
        'STREAMobj must be a single stream, i.e. there is only one channelhead.');
end
