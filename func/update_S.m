
function S = update_S(S,D,X,U)

H = U*U';
G = D'*D;
L = D'*X*U';


L_plus  = (abs(L)+L)/2;
L_minus = (abs(L)-L)/2;

G_plus  = (abs(G)+G)/2;
G_minus = (abs(G)-G)/2;

H_plus  = (abs(H)+H)/2;
H_minus = (abs(H)-H)/2;

term_1 = H_minus*S*G_plus + H_plus*S*G_minus  + L_plus  + eps;
term_2 = H_plus*S *G_plus + H_minus*S*G_minus + L_minus + eps;

S = S.*sqrt(term_1./term_2);
S = S*diag(sqrt(1./(diag(S'*S)+eps)));