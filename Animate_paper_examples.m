
clear all
close all
clc
warning off

%% make names
dirname = 'Examples_in_paper\';

name{1}	= [dirname,'Single_link_1']; %ok
name{2} = [dirname,'Single_link_2']; %ok
name{3} = [dirname,'Double_link_1']; %ok
name{4} = [dirname,'Double_link_2']; %ok
name{5} = [dirname,'Double_link_transmission_1']; %ok
name{6} = [dirname,'Double_link_transmission_2']; %ok
name{7} = [dirname,'Triple_link_1']; %ok
name{8} = [dirname,'Triple_link_2']; %ok

%% animate all examples from paper and save .gif in Animations folder
for k =  1:8
    close all
    load([name{k},'.mat'])
    animateDNA(DNA,t,state)
end

