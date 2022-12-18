

function A =  get_ram_dicionary_indicator_mtrx(periods)

% e.g A = [1,0 0 0 0 0 ;
%      0,1,0 0 0 0
%      0,0,1,1,0,0]

P_max =  max(periods);
LEN   = length(periods);

A =  zeros(P_max,LEN);

for i =  1:P_max
    index_tmp =  find(periods==i);
    A(i,index_tmp) = 1;
end

