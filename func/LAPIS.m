
function  [Factors,detected_periods,running_time] = LAPIS(Y,opts)

% INPUT:
% Y: observation with missing values
% opts:

% OUTPUT:
% X: completed observation


% for period dictionary
if ~isfield(opts, 'Dictionary_type'),   opts.Dictionary_type  = 'Ramanujan'; end
if ~isfield(opts, 'Pmax'),              opts.Pmax             = 50;          end
if ~isfield(opts, 'A'),                 opts.D                = Create_Dictionary(opts.Pmax,size(Y,1),opts.Dictionary_type); end
if ~isfield(opts, 'Penalty_type'),      opts.Penalty_type     = 'square';    end


if ~isfield(opts, 'halting_Thr'),       opts.halting_Thr = 1e-3; end  % for algo
if ~isfield(opts, 'max_iter'),          opts.max_iter    = 100;  end  % for algo

if ~isfield(opts, 'visual'),            opts.visual      = 0;   end  %

if ~isfield(opts, 'rho'),               opts.rho         = 1.1; end
if ~isfield(opts, 'is_show'),           opts.is_show     = 0;   end

if ~isfield(opts, 'lambda_1'),          opts.lambda_1    = 0.1; end
if ~isfield(opts, 'lambda_2'),          opts.lambda_2    = 0.1; end
if ~isfield(opts, 'lambda_3'),          opts.lambda_3    = 0.5; end


Pmax        = opts.Pmax;
halting_Thr = opts.halting_Thr;
max_iter    = opts.max_iter;

%% get penalty matrix and indicator matrix
[periods,penalty_vector]  =  get_period_penalty(Pmax,opts.Penalty_type);
H_inv                     =  diag(1./penalty_vector);                   % inverse of diagnal matrix
A                         =  get_ram_dicionary_indicator_mtrx(periods); % capture block structure


QUIET   = 0;
t_start = tic;


%% initialization
D      =  opts.D* H_inv;      % H is diagnal penlty matrix
W      =  double(logical(Y)); % mask of x: 1/0, binary vector

[m,n]  = size(D);


S = eye(n,n);
X = Y;


disp('before loop')
%% main loop from here
iter = 0;
while iter < max_iter

    iter =  iter+1;
    disp(iter)

    D_hat     = D*S;
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %++++++++++++++++++++++ update U
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    [U,U_plus,U_minus] = update_MPID_U(X,A,D_hat,opts);

    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %++++++++++++++++++++++ update X
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    X = update_X(W,Y,D_hat,U,opts);

    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %++++++++++++++++++++++ update S
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    S = update_S(S,D,X,U);


    %++++++++++++++++++++++ check convergence
    if iter > 1
        history.objval(iter)  = getObj(X,Y,W,D_hat,A,U,U_plus,U_minus,opts);
        obj_resi              = abs(history.objval(iter)- history.objval(iter-1))/history.objval(iter-1);

        if mod(iter,5) == 0
            disp(['Out loop: Iter:',num2str(iter),'::','obj_res = ',num2str(abs(obj_resi))]);
        end

        if obj_resi <  halting_Thr
            break;
        end
    end

end

%% output

Factors.X = X;
Factors.S = S;
Factors.U = U;
Factors.U_plus = U_plus;
Factors.U_minus = U_minus;


if ~QUIET
    running_time = toc(t_start);
end


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% for EACH sample period U
detected_periods = zeros(size(U,2),Pmax);
for i = 1:size(U,2)
    Tmp_period            =  U(:,i);
    s_individual          =  H_inv * Tmp_period;
    periods_vector_tmp    =  compute_vector_period(s_individual,opts);

    [periods_output_mag, periods_idx_tmp] = sort(periods_vector_tmp,'descend');
    detected_periods(i,:)                 = periods_idx_tmp;% (1:TOP_K); % top-1

    period_vector_all(i,:) =  periods_vector_tmp;
end

Factors.period_vector_all = period_vector_all;

end



function obj = getObj(X,Y,W,D_hat,A,U,U_plus,U_minus,opts)
obj = 0.5* norm(X-D_hat*U,'fro') + opts.lambda_1 * norm(U,1) ...
    + opts.lambda_2 * nuclear_norm(A*(U_plus+U_minus)) + opts.lambda_3*norm(W.*(X-Y),1);
end




function detect_periods_vector = compute_vector_period(s,opts)


Pmax              = opts.Pmax;
energy_s          = 0.*[1:Pmax];
current_index_end = 0;

for i=1:Pmax
    k_orig = 1:i;k=k_orig(gcd(k_orig,i)==1);
    current_index_start = current_index_end + 1;
    current_index_end   = current_index_end + size(k,2);

    for j = current_index_start:current_index_end
        energy_s(i) = energy_s(i)+((abs(s(j)))^2);
    end
end

energy_s(1)    = 0;
detect_periods_vector =  energy_s./norm(energy_s);
end

