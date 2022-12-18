


function [data_matrix,statis] = build_syn_data_matrix(opts_multivariate)

% every time series has shared  period
% randomly select a few time series, give each of them a unique period
if ~isfield(opts_multivariate, 'incomplete'),          opts_multivariate.incomplete          = 0; end  %


num_groups       =  opts_multivariate.num_groups;
num_smaple       =  opts_multivariate.num_smaple;
Input_Datalength =  opts_multivariate.Input_Datalength;
syn_data_matrix  =  zeros(Input_Datalength,num_smaple);


% randomly sort
% individual_period_all      = opts_multivariate.Input_individual_periods{1};
% individual_period_all_rand = individual_period_all(randperm(length(individual_period_all)));


%% settings for univariate

SNR                   =  opts_multivariate.SNR;
% Input_periods  =  opts_multivariate.Input_periods{1};


% for i = 1:num_smaple
%     individual_period_tmp  =  individual_period_all_rand(i);
%     period_tmp             =  [Input_shared_periods,individual_period_tmp];
%     syn_data_matrix(:,i)   =  generate_periodic_signal(period_tmp,Input_Datalength,SNR);
%
%     statis.period{i}       =  period_tmp; % period for each time series
% end
groupIndex        = get_group_index(num_smaple,num_groups);
statis.groupIndex = groupIndex;

for nn = 1:num_groups
    IDX        = find(groupIndex==nn);
    period_tmp = opts_multivariate.Input_periods{nn};

    statis.group_period{nn} = period_tmp;

    for i = 1:length(IDX)
        id_tmp = IDX(i);
        syn_data_matrix(:,id_tmp)   =  generate_periodic_signal(period_tmp,Input_Datalength,SNR);
        all_period{id_tmp,1}        =  period_tmp;  % period for each time series
        all_period{id_tmp,2}        =  nn; % group ID
    end
end

statis.period =  all_period;

if opts_multivariate.visual_signal
    for i = 1:num_smaple
        x = syn_data_matrix(:,i) ;
        figure, plot(x);
        title('Input Signal Without missing data','FontSize',18);
        xlabel('Time index','FontSize',18);
        ylabel('Signal amplitude','FontSize', 18)
        set(gca,'looseInset',[0 0 0 0])
    end
end



statis.raw_complete_data_matrix = syn_data_matrix;
%% random drop values in data matrix
if opts_multivariate.incomplete

    ratio_incomplete   = opts_multivariate.ratio_incomplete;
    num_total_elements = num_smaple * Input_Datalength;
    num_missing        = floor(num_total_elements*ratio_incomplete);

    missing_idx        = randperm(num_total_elements,num_missing);
    data_vector        = syn_data_matrix(:);

    data_vector(missing_idx) = 0;
    data_matrix              = reshape(data_vector,[Input_Datalength,num_smaple]);

    for i = 1:num_smaple
        statis.missing_ratio{i}  = nnz(~data_matrix(:,i)) /Input_Datalength;
    end


    if opts_multivariate.visual_incomplete
        for i = 1:num_smaple
            x_incomplete = data_matrix(:,i) ;
            figure, plot(x_incomplete);
            title('Input Signal With Missing data','FontSize', 18);
            xlabel('Time index','FontSize', 18);
            ylabel('Signal amplitude','FontSize', 18)
            set(gca,'looseInset',[0 0 0 0])
            %  set(h.a2, 'ButtonDownFcn',@buttonDownCallback)
        end
    end
else

    data_matrix=  syn_data_matrix;
end



end


% function M = assign_period_group(opts_multivariate,num_groups,groupIndex,num_period_each)
%
% Input_periods = opts_multivariate.Input_periods{1};
% totoal_period = length(Input_periods);
% % periods_idx   = ranperm(totoal_period);
% % Input_periods = Input_periods(periods_idx);
%
% period_matrix = zeros(num_groups,num_period_each);
%
% for i = 1:num_groups
%     id_tmp = randample(totoal_period,num_period_each);
%
%     M = reshape(Input_periods,[num_groups,totoal_period/num_groups]);
%
%     id_period = [id_period,id_tmp]
% end
%
% end






