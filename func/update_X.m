
function   X = update_X(W,Y,D_hat,U,opts)



% for ADMM
%if ~isfield(opts, 'eta'),                 opts.epsilon =1e-4; end  %
if ~isfield(opts, 'tau'),                 opts.tau     =1.1; end  %
if ~isfield(opts, 'rho'),                 opts.rho     =1.1; end  %
if ~isfield(opts, 'MAX_mu'),              opts.MAX_mu =1e10; end  %
if ~isfield(opts, 'max_max_iter_beta'),   opts.max_max_iter_beta = 100; end  %


% epsilon = opts.epsilon ;
% tau     = opts.tau;
% max_rho = opts.max_rho ;

% max_iter = opts.max_max_iter_beta;

QUIET    = 0;
MAX_ITER = 50;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;
MAX_mu   = 10^10;
mu       = 1; %10^-3;
rho      = 1.1;

% I_L = eye(size(L));
% Temp_1 =  D'*D-lambda_1*L-mu*I_L;

% if ~QUIET
%     fprintf('%3s\t%10s\t%10s\n', 'iter', 'X norm', 'E pri', 'objective');
% end


DU   =  D_hat*U;
%   I_G =  eye(size(G));
%Temp_1 = (I_G+2*opts.lambda_3*G+ 0.5*mu*I_G)^(-1);

E =  zeros(size(Y));
M =  zeros(size(Y));
%%
for k = 1: MAX_ITER
    
    
    %++++++++++++++++++++ U-update
    H     = E + Y + M/mu;
    X     = (DU+ mu*H)/(1+mu);
    
    %++++++++++++++++++++ V-update
    T     = X-Y- M/mu;
    E     = max(abs(T)-(opts.lambda_3/mu).*W,0).*sign(T);
    
    %++++++++++++++++++++ M-update
    M = M + mu*(E- (X - Y));
    
    %++++++++++++++++++++ mu-update
    % mu = max(mu * rho,MAX_mu);
    
    % check convergence
    if k >1
        history.objval(k)  = get_obj(X,Y,W,DU, opts);
        
        history.X_norm(k)  = norm(X-X_old);
        history.E_norm(k)  = norm(E-E_old);
        
%         if ~QUIET
%             fprintf('%3d\t%10.2f%10.2f%10.2f\n', k, history.X_norm(k), history.E_norm(k),  history.objval(k));
%         end
        if (history.X_norm(k) < ABSTOL && history.E_norm(k)<ABSTOL && history.objval(k) < ABSTOL)
            break;
        end
        
    end
    X_old =  X;
    E_old =  E;
    
end

end


function obj = get_obj(X,Y,W,DU, opts)
obj = 0.5*norm(X-DU,'fro') + opts.lambda_3*norm(W.*(X-Y),1);
end




