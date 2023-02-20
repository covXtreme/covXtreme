function beta = ipqr (q, y, X)
%   Use interior point algorithm for quantile regression.
%   (See Koenker and Park (1996), X. Econometrics) 
%
%  Inputs:
%     -- q is the quantile desired (number in (0,1)).
%     -- y is the response vector (n X 1).
%     -- X is the predictor matrix (n X m).

eta=.97;   % Value recommended by K&P
toler=1e-5;

p=size(X,2);
beta=quantile(y,q)*ones(p,1);
iteration = 0;
n=length(y);
d=zeros(n,1);
change = realmax;

L=1e-2;  %ridge penalty

residuals = y-X*beta;      
lastobj=q*sum(residuals)-sum(residuals(residuals<0));

while abs(change) > toler 
	iteration = iteration + 1;    
	dhat=d-X*((X'*X+L*eye(p))\(X'*d));
	m=toler+max([-dhat./(1-q);dhat./q]);
	d = dhat ./ m;

    % Now we perform the inner iterations.  K&P recommend two of them per
    % outer iteration.
    for k=1:2
        dmatrix = min([q-d  1-q+d],[],2);
        d2=dmatrix.^2;
        direc =(X'*(d2(:,ones(p,1)).*X)+L*eye(p))\(X'*(d2.*residuals));
        s = d2.*(residuals - X*direc);
        alpha = max(max([s./(q-d) -s./(1-q+d)]));
        d=d+eta/alpha*s;
    end
    
    % Now we do the line search recommended by K&P.  Some of the parameters
    % passed to the fmin function are ad hoc (K&P give no specific
    % recommendations).
    step=fminbnd(@(t)ipqr_objfunc(t,q,y,X,beta,direc),-1,1,[]) * direc;
    
    beta = beta + step;
    residuals = y - X*beta;
    obj=q*sum(residuals)-sum(residuals(residuals<0));
    change=lastobj-obj;
    lastobj=obj;
end

function value = ipqr_objfunc (t,q,y,X,start,dir)
% Calculate objective function at start+t*dir. 
res=y-X*(start+t*dir);  % Find residuals for proposed beta.
value=q*sum(res) - sum(res(res<0));
