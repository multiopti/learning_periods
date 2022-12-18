

function [x,x_incomplete] =  data_syn_generator(opts_syndata)

rng('default')

% opts_syndata.Input_Datalength = 200;


%%


% Step 0: Selection of Dictionary Parameters
%
% Pmax = 50; %The largest period spanned by the NPDs
% Dictionary_type = 'Ramanujan'; %Type of the dictionary
%


% --- signal parameters ++++++++++++++
if ~isfield(opts_syndata, 'Input_Datalength'),  opts_syndata.Input_Datalength   = 200; end  %
if ~isfield(opts_syndata, 'Input_Periods'),     opts_syndata.Input_Periods      = {[3,7,11]}; end
if ~isfield(opts_syndata, 'SNR'),               opts_syndata.SNR                = 10; end  % Signal to noise ratio of the input in dB
if ~isfield(opts_syndata, 'visual_incomplete'), opts_syndata.visual_incomplete  = 0; end
if ~isfield(opts_syndata, 'visual_signal'),     opts_syndata.visual_signal      = 0; end


x   = generate_periodic_signal(opts_syndata.Input_Periods{1},opts_syndata.Input_Datalength,opts_syndata.SNR);

if opts_syndata.visual_signal
    figure, plot(x);
    title('Input Signal Without missing data', 18);
    xlabel('Time index','FontSize', 18);
    ylabel('Signal amplitude','FontSize', 18)
    set(gca,'looseInset',[0 0 0 0])
    
end

%---- incomplete signal ++++++++++++++
if ~isfield(opts_syndata, 'incomplete'),          opts_syndata.incomplete          = 0; end  %
if ~isfield(opts_syndata, 'ratio_incomplete'),    opts_syndata.ratio_incomplete    = 0.2; end  %
if ~isfield(opts_syndata, 'missing_window_size'), opts_syndata.missing_window_size = 1; end  % Signal to noise ratio of the input in dB
if ~isfield(opts_syndata, 'visual_incomplete'),   opts_syndata.visual_incomplete   = 0; end  % Signal to noise ratio of the input in dB

if opts_syndata.incomplete
    x_incomplete = gen_incomplete(x,opts_syndata.Input_Datalength ,opts_syndata.ratio_incomplete, opts_syndata.missing_window_size);
    x_incomplete = x_incomplete./norm(x_incomplete,2);
    if opts_syndata.visual_incomplete
        figure, plot(x_incomplete);
        title('Input Signal With Missing data','FontSize', 18);
        xlabel('Time index','FontSize', 18);
        ylabel('Signal amplitude','FontSize', 18)
        set(gca,'looseInset',[0 0 0 0])
        %  set(h.a2, 'ButtonDownFcn',@buttonDownCallback)
    end
else
    
    x_incomplete = [];
end







