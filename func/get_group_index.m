

function groupIndex = get_group_index(NSample,NG)
% NSample: total sample
% NG: groups


% Your matrix
% NSample = 50;
% NCOL = 100;

% M = rand(NSample,NCOL);
% NG = 6; % The number of groups

% Generate an index that assigns each row to one of the groups
groupIndex = randi(NG,NSample,1);


% % Assign the rows to their own matrix (in a cell array)
% Mg = cell(NG,1);
% for nn = 1:NG
%     Mg{nn} = M(groupIndex==nn,:);
% end
