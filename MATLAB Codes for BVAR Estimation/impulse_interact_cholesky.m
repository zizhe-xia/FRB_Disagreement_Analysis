function response=impulse_interact_cholesky(By,Bx,Bw,smat,nstep,x,full,exo_included)
% function response=impulse_interact(By,Bx,smat,nstep,x), C. Sims' code. 
% Modified by Zizhe, Xia, 2018
% The value of x is the pre-set exogenous variable to calculate the
% interaction term.
% The By matrix is the relevant coefficent matrix, i.e. part of the coefficient
% matrix from a single draw which contains the variables considered in the
% impulse response function; Bx is the drawn coefficient for the
% interaction terms. 
% By is a neq x nvar x nlags matrix.  neq=nvar, of course, but the first index runs over 
% equations. In response, the first index runs over variables, the second over 
% shocks (in effect, equations).
% IF full=true /  1, then interaction is allowed for x and all Y
% IF full=false / 0, then interaction is allowed only for x and the last
% variable which is effective FFR
if full
    [neq,nvar,nlag]=size(By);
    response=zeros(nvar,neq,nstep);
    response(:,:,1)=smat'; % need lower triangular, last innovation untransformed
    for it=2:nstep
       for ilag=1:min(nlag,it-1)
           response(:,:,it)=response(:,:,it)+(By(:,:,ilag) + x*Bx(:,:,ilag))*response(:,:,it-ilag);
           % By(:,:,ilag)*response(:,:,it-ilag)+x*Bx(:,:,ilag)*response(:,:,it-ilag);
       end
    end
else
    [neq,nvar,nlag]=size(By);
    response=zeros(nvar,neq,nstep);
    response(:,:,1)=smat'; % need lower triangular, last innovation untransformed
    for it=2:nstep
       for ilag=1:min(nlag,it-1)
          response(:,:,it)=response(:,:,it)+(By(:,:,ilag) + [x*Bx(:,:,ilag) zeros(nvar,nvar-1)])*response(:,:,it-ilag);
          % By(:,:,ilag)*response(:,:,it-ilag)+x*Bx(:,:,ilag)*response(:,:,it-ilag);
       end
    end
end

% Original response function using decomposed errors as the shock
function response=impulse(By,smat,nstep)
% function response=impulse(By,smat,nstep), C. Sims' code.
% smat is a square matrix of initial shock vectors.  To produce "orthogonalized
% impulse responses" it should have the property that smat'*smat=sigma, where sigma
% is the Var(u(t)) matrix and u(t) is the residual vector.  One way to get such a smat
% is to set smat=chol(sigma).  To get the smat corresponding to a different ordering,
% use smat=chol(P*Sigma*P')*P, where P is a permutation matrix.
% By is a neq x nvar x nlags matrix.  neq=nvar, of course, but the first index runs over 
% equations. In response, the first index runs over variables, the second over 
% shocks (in effect, equations).
[neq,nvar,nlag]=size(By);
response=zeros(nvar,neq,nstep);
response(:,:,1)=smat'; % need lower triangular, last innovation untransformed
for it=2:nstep
   for ilag=1:min(nlag,it-1)
      response(:,:,it)=response(:,:,it)+By(:,:,ilag+1)*response(:,:,it-ilag);
   end
end
