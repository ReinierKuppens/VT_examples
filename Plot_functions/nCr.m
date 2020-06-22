%   C = nCr(n,2) returns the number of combinations with
%   
%   n the number of values to pick from 
%   k the number of picks
%
%   with the formula: n!/(k!(n-k)!)
%
%   Be aware that the factorial function can only be used up to n = 170, 
%   because of machine precision. 170!=7.2574*10^306 to give an idea 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   Written by Reinier Kuppens as part of Msc Thesis                      %           
%   Original Date: 26-05-2015                                             %
%   Last modified: 29-06-2015                                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = nCr(n,k)

if n>100
    error('nCr contains a list with 100 elements. Please select n<100 or add more elements to the list')
end

%See https://oeis.org/A000217 for the sequence.

A000217=[  0
           1
           3
           6
          10
          15
          21
          28
          36
          45
          55
          66
          78
          91
         105
         120
         136
         153
         171
         190
         210
         231
         253
         276
         300
         325
         351
         378
         406
         435
         465
         496
         528
         561
         595
         630
         666
         703
         741
         780
         820
         861
         903
         946
         990
        1035
        1081
        1128
        1176
        1225
        1275
        1326
        1378
        1431
        1485
        1540
        1596
        1653
        1711
        1770
        1830
        1891
        1953
        2016
        2080
        2145
        2211
        2278
        2346
        2415
        2485
        2556
        2628
        2701
        2775
        2850
        2926
        3003
        3081
        3160
        3240
        3321
        3403
        3486
        3570
        3655
        3741
        3828
        3916
        4005
        4095
        4186
        4278
        4371
        4465
        4560
        4656
        4753
        4851
        4950];

C = A000217(n);
    
end