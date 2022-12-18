
function  [U,U_plus,U_minus] = update_MPID_U(X,A,D_hat,opts)

% min 0.5||X-DSU||_F^2 + lambda_1 (||A*U_plus||_1+||A*U_minus||_1) + lambda_2 ||A(U_plus+U_minus)||_*
%

lambda_1 =  opts.lambda_1;
lambda_2 =  opts.lambda_2;

max_iter    = 50;
halting_Thr = 1e-4;



%% initialization
[p,L] = size(A);
[T,N] = size(X);


P_1 = zeros(L,N);
P_2 = zeros(L,N);


Theta_1 = zeros(L,N);
Theta_2 = zeros(p,N);


rho_1  = 1e-6;
rho_2  = 1e-6;

tau    = 1.1;

% D_hat =  D*S;
DtX   =  D_hat'*X;
DtD   =  D_hat'*D_hat;
AtA   =  A'*A;

I     =  eye(L,L);

%% main loop

for iter =  1:max_iter
    
    %++++++++++++++++++++++ update  U_plus
    
    temp_1 = DtD + rho_1* I;
    temp_2 = DtX + rho_1*P_1 + Theta_1;
    % U      = temp_2\temp_1;
    U      = inv(temp_1)*temp_2;
    
    U_plus  = (abs(U)+U)/2;
    U_minus = (abs(U)-U)/2;
    
    %++++++++++++++++++++++ update proxy: P1,,P2, P3
    F_1 = U - Theta_1/rho_1;
    P_1 = max(abs(F_1)-lambda_1/rho_1,0).*sign(F_1);
    
    
    F_2 = A*(U_plus + U_plus) - Theta_2/rho_2;
    [L,Sigma,R] = svd(F_2,'econ');
    Sigma_thr   = diag(max(diag(Sigma) - lambda_2/rho_2 ,0));
    P_2         = L*Sigma_thr*R';
    
    %++++++++++++++++++++++ update  Lagrange multipliers: Theta
    Theta_1 = Theta_1 + rho_1 *(P_1 - U);
    Theta_2 = Theta_2 + rho_2 *(P_2 - A*(U_plus + U_plus));
    
    
    % ++++++++++++++++++++++ update penalty parameters
    rho_1  = min(tau*rho_1,1e10);
    rho_2  = min(tau*rho_2,1e10);
    

    
    %++++++++++++++++++++++ check convergence
    if iter > 1
        % iter
        history.objval(iter)  =  getObj(X,D_hat,A,U,U_plus,U_minus,opts);
        obj_resi =  abs(history.objval(iter)- history.objval(iter-1))/history.objval(iter-1);
        
        if mod(iter,5) == 0
            disp(['Iter:',num2str(iter),'::','obj_res = ',num2str(abs(obj_resi))]);
        end
        
        if obj_resi <  halting_Thr
            break;
        end
    end
    
    
    if iter == max_iter
        fprintf('UPDATE U: MAximum number of iteration (%g) reached! \n', iter)
    end
end


end



function  obj = getObj(X,D_hat,A,U,U_plus,U_minus,opts)

obj = 0.5* norm(X -D_hat*U,'fro')^2 + opts.lambda_1 * norm(A*U,1) ...
    + opts.lambda_2 * nuclear_norm(A*(U_plus+U_minus));

end
