function ...                                % RSLEASQR.M   P. 1 of 4
   [ pbest,se,dep,  rbest,fbest,sbest,  iter,cnvrg,qi ] = ...
      rsleasqr( resid,drdp,  p,dpr,dpa,ponoff, ...
      stol,niter,erra,  varargin );
% rsleasqr f: Fit a nonlinear function to data.
%         [ p, se, dep, ...             % Parameters & error
%           r, f, ss, ...               % Fitted function. 
%           iter, cnvrg, qi ] = ...     % Convergence & covariance.
%           leasqr( resid, drdp, ...    % Function names.
%           pin, dpr, dpa, ponoff, ...  % Parameters & increments,
%           stol, niter, erra, ...      % Convergence & error control.
%           P1, P2, ... );              % User variables.
% See header comments in listing for details.

% ------------INPUT VARIABLES (all vectors are columns) --------------
% resid = string name of the function .m file.  See example.
%        The function in resid should have the calling sequence:
%           [ ss, r, f, B1, B2, ... ] = resid( p, A1, A2, ... ), where
%           p  = the parameters to be adjusted.
%           A1, A2 ... = miscellaneous user's input arguments
%                        such as x, y, w, et al.
%           f  = the computed values.
%           r  = w .* (y-f),    where   y = data, and
%                w = square roots of weights (if any).
%           ss = r(:)'*r(:), the sum of squares.
%           B1, B2 ... = miscellaneous user's output arguments.
%         The function resid may abort the run by setting ss<0.
%         ss = -1 signals a "soft" abort, meaning that
%          rbest and fbest are correct, but some user-defined bound
%          has been exceeded, usually by pbest or fbest, so that 
%          the estimate should not be accepted as a solution.
%          However, rbest can be used in internal calculations.
%         ss = -2 signals a "hard" abort, meaning that
%          fbest and rbest could not be computed, e.g. beacuse 
%          logs or sqrts of negative numbers would have occured.
%         Despite these signals, leasqr will make some effort to
%         find p that produces proper results, but if this proves
%         impossible, e.g. if initial p produces an abort, 
%         then output ss (see sbest below) will be -1 or -2.
% drdp = string name of the partials .m file.  See example.
%        Computes the matrix jac of partial derivatives of the residual
%        vector described above,  i.e. jac(i,j) = dr(i)/dp(j).
% p = vector of initial parameters to be adjusted by leasqr.
% dpr = relative increments of p for numerical partials (see also dpa).
%       I.e. the increment used will be abs(dpr(j)*p(j)).
%       dpr(j) > 0 means use central differences.
%       dpr(j) = 0 means use absolute increment dpa(j).   BUT IF dpa(j)
%                  IS ALSO 0, HOLD p(j) FIXED, SET ITS PARTIALS TO 0,
%                  EVEN IF NUMERICAL PARTIALS ARE NOT USED!!
%       dpr(j) < 0 means use one-sided differences.
% dpa = absolute increments of p for numerical partials.   Often, this
%       is a guard against minuscule relative increments when p(j)
%       becomes small.   If both dpr(j) and dpa(j) are nonzero, the
%       increment used is max { abs[dpr(j)*p(j)], abs[dpa(j)] }.
%       If dpr(j)~=0, the sign of dpr(j) governs the use of central
%       or one-sided differences.   If dpr(j)=0, the sign of dpa(j)
%       governs that choice in the same manner.   If dpa(j)=0, and
%       dpr(j)~=0, relative increments are always used.   If both dpr(j)
%       and dpa(j) = 0, p(j) is held fixed, its partials set to 0.
% ponoff = vector of p activity. If ponoff(i)=0, p(i) is held fixed.
%          If ponoff is empty, all p are active (not fixed).
% stol = scalar tolerance on fractional improvement in the sum of squares.
%        let ss = current sum of squares, and ssprev = ss from previous
%        iteration.  Then convergence is declared if sprev-ss<ss*stol.
% niter = scalar maximum number of iterations, i.e. calls to drdp.
%         if niter = 0, outputs are given for the input values of p.
% erra = 1 for error analysis option, 0 otherwise.
% varargin = a set of optional arguments which, if provided, are passed to
%            the functions resid and drdp.   This option enables the user
%            to communicate with his functions without global variables.
% ------------OUTPUT VARIABLES----------------- % LEASQR.M   P. 3 of 4
% pbest = vector of final parameters, i.e. the solution.
% se    = asymptotic standard errors of p.
% dep   = dependency factors for p; dep(j) = ratio of se(j) to that se(j)
%         which would have been produced if p(j) were the only parameter.
% rbest = matrix of best residual values computed by resid.
% fbest = matrix of best function values computed by resid.
% sbest = scalar sum of squares = sum-over-i(wt(i)*(y(i)-f(i)))^2.
%         if sbest<0, the function resid has aborted the run.
%         See input function resid about soft & hard aborts.
% iter  = scalar number of iterations used.
% cnvrg = scalar = 1 if convergence, = 0 otherwise.
% qi    = matrix covariance of p.
% ----------- EXAMPLE. Begin hyperb --------------------------------------
%   function [ss,r,f] = hyperb(p,x,y)  % hyperb = function to be fitted
%   vmax = p(1);   km = p(2);          %   by adjusting p = [ vmax, km ].
%   f = vmax*x./(km+x);                % f = fitting function.
%   r = y-f;                           % r = residuals.
%   ss = sum(r.*r);                    % ss = sum of squares.
% ----------- End hyperb. Begin dhyperb ----------------------------------
%   function dr = ...
%      dhyperb(hyperb,p,dpr,dpa,r,x,y)          % dhyperb computes
%   vmax=p(1);   km=p(2);                       %  dr(i)/dp(j), i.e.
%   dr = [ -x./(km+x), vmax*x./((km+x).^2) ];   %  dr/dvmax and dr/dkm.
% ----------- End dhyperb. Begin master session --------------------------
%   vmax=1;  km=1;  p=[vmax;km];             % Simulation study. First set
%   x=(0:.1:2)';    y=zeros(21,1);           % p to known values, then
%   [s,r,f]=hyperb(p,x,y);                   % generate artificial data.
%   y=f+.01*randn(21,1);                     % Add noise to create "data".
%   p = [ .9*vmax; 1.1*km ];  niter=50;      % Perturb p.  See how well p
%   dpr=ones(2,1); dpa=dpr; stol=.001;       % is recovered by fitting to
%   [r,p,se,dep,r,f,ss,iter,cnvrg]= ...      % the noisy data.  Note: x & y
%   leasqr('hyperb','dhyperb', ...           % are arguments to be passed
%      p,dpr,dpa,stol,niter,x,y);            % to hyperb & dhyperb.
%   disp('   Parameters    Std. errs.  Dependencies');
%   disp([p,se,dep]);                        % Show parameters et al.
%   plot(x,y,'o',x,f,'-');                   % Plot y and fitted f.
% ---------- End master session. End EXAMPLE. ---------------------------
   % ------------BEGIN PROGRAM---------------- % RSLEASQR.M   P.4 of 4
   ID='rsleasqr';
   np=length(p);   p=reshape(p,np,1);         % How many parameters?
   pbest=p;                                   % Initalize best p.
   se=[];   dep=[];   qi=[];                  % Empty se, dep, qi.
   iter = 0;   cnvrg = 0;                     % Initial iter, cnvrg.
   dpr=reshape(dpr,np,1);                     % Relative increments.
   dpa=reshape(dpa,np,1);                     % Absolute increments.
   zcol=zeros(np,1);   zrow=zeros(1,np);      % Vectors for fixing p.
   if isempty(ponoff);                        % If ponoff is empty,
      ponoff = ones(size(p));   end;          %  make all p active.
   fixp=((dpr==0)&(dpa==0))|(~ponoff);        % Any p held fixed?
   sc = zeros(np,np);   id = eye(np);         % Scaling matrices.
   epstab=[ 1e-10, .001, .1, 1, 10, 100 ];    % Table of trial epsilons.
   [sbest, rbest, fbest] = ...
      feval( resid, p, varargin{:} );         % Initialize best ss, r, f.
   if sbest<0; return; end;                   % Abort if sum sqs. <0.
   [nrr,ncr]=size(rbest);  nr=nrr*ncr;        % Dimensions of r.
   for iter=1:max([1,niter]),  pprev=pbest;   % Iteration loop.
     prt = feval( ...                         % Partial derivatives.
       drdp, resid, p, dpr, dpa, ...
       ponoff, rbest, varargin{:} );
     q=prt'*prt;                              % Inner product.
     sprev=sbest;  sgoal=(1-stol)*sprev;      % Save best & goal ss.
     r=rbest;   g=prt'*reshape(r,nr,1);       % Best residuals.
     for j=1:np,                              % For each p:
       if fixp(j) | q(j,j)==0,                %   Scan for fixed p.
         q(:,j)=zcol;   q(j,:)=zrow;          %   Zero out row & col. &
         q(j,j)=1;   g(j)=0;                  %     insure nonsingular q,
       else   sc(j,j)=1/sqrt(q(j,j));         %     or scale q to avoid
       end;  end;                             %     "bad conditioning"
     q=sc*q*sc;   g=sc*g;                     %     messages.
     if niter==0, iter=0; break; end;         % Error analysis only.
     for epsl=epstab,                         % Trial epsilon loop.
       p=pprev-sc*((q+epsl*id)\g);            % Trial parameters.
       [ss,r,f] = ... 
         feval(resid,p,varargin{:});          % Sum sqs, resids, & func.
       if ss>=0;                              % Trial not aborted.
         if ss<sbest,  sbest=ss;              % If any improvement,
           pbest=p;  rbest=r;                 %   save best values of
           fbest=f;  end;                     %   ss, p, r, and f.
         if ss<=sgoal, break; end;            % If small ss, stop eps loop
       end; end;
     if sbest>sgoal | sbest==0,               % If sbest is big or 0,
       cnvrg=1;   break;   end;   end;        %   we have converged.
   if erra==0;   return;   end;               % No error analysis.
   qi=inv(q+1e-14*id);                        % Error analysis.
   sef=sqrt(sbest/ ...
     (max([1,nr-np+sum(fixp)])));             % Std. error of the fit.
   ze=~fixp;                                  % Which p are not fixed?
   if max(size(q))==1,                        % Std.errors,
      se=sef*sc*sqrt(qi)*ze;   dep=ze;        %   dependencies, and
   else se=sef*diag(sc).*sqrt(diag(qi).*ze);  %   covariance of
     dep=sqrt(diag(qi).*ze.*diag(q));   end;  %   active p.
   qi = qi*sef^2;                             % Scale qi.
   for j=1:np;                                % Zero out qi
      if fixp(j);   qi(j,j)=0;   end;   end;  %  for fixed p.
