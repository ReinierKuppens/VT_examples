% plotGrounds(DNA) plots all the grounds of a mechanism.
%
%   Input:  DNA
%   Output: figure with all grounds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 2-07-2015                                              %
%   Last modified: 2-07-2015                                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function []=plotGrounds(DNA)
hold on

% cdata       = [1 1 1];

[xOutline,yOutline,xSmallCircle,ySmallCircle,xlines,ylines] = getGroundOutline2(DNA);

wLine = 1;


for k = 1:numel(xOutline)
    fill(xOutline{k},yOutline{k},[1 1 1],'Linewidth',wLine);
    fill(xSmallCircle{k},ySmallCircle{k},'w','Linewidth',wLine);
    plot(xlines{k},ylines{k},'k','linewidth',wLine);

end










