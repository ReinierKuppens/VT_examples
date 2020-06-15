%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %
%   Original Date: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function DNA=initializeDNA

DNA.incstr = [];
DNA.edgelabel = [];
DNA.Mpar  = zeros(0,3);
DNA.Hpar  = zeros(5,0);
DNA.Spar  = zeros(6,0);
DNA.Ppar  = zeros(3,0);
DNA.Pmpar = zeros(5,0);
DNA.fitness = 0;
DNA.bestTrajectory =zeros(0,2);