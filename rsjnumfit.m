function prt=rsjnumfit( ...                           rsjnumfit   p.1/2
   resid, p, dpr, dpa, ponoff, r, varargin );
% rsjnumfit: Numerical partial derivatives dr/dp for use with leasqr.m;
%         PLEASE NOTE: ALL INPUT VECTORS ARE COLUMNS!!!
% ------------ INPUT -----------------------------------------------------
% resid = string naming the residual .m file, which computes the
%         residual matrix r = data minus computed function or the
%         negative, optionally weighted, i.e. w.*(y-f) or w.*(f-y).
%         Only the 1st 2 output arguments of resid, namely ss and r,
%         are used here.   If ss=-1 (soft abort), the calculation 
%         proceeds, but if ss=-2 (hard abort), the resulting r cannot
%         be used.   Some attempt will be made to avoid the problem
%         by using 1-sided instead of 2-sided differences, or by
%         reversing the sense of 1-sided differences.   If none of
%         these strategies work, the offending column of prt (see below)
%         is set to zero.   
% p     = column of current parameter values.
% dpr   = column of relative increments of p for numerical partials.
%         dpr(j)>0 means use 2-sided differences.
%         dpr(j)=0 means use absolute increment dpa(j).   But if dpa(j)
%                  is also 0, hold p(j) fixed, set its partials to 0.
%         dpr(j)<0 means use 1-sided differences.
% dpa   = column of absolute increments of p for numerical partials.
%         Often, these guard against an increment getting too small when
%         p(j) becomes small.   If both dpr(j) and dpa(j) are nonzero,
%         the increment used is max{ abs[dpr(j)*p(j)], abs[dpa(j)] }.
%         If dpr(j)=0, 
%         dpa(j)>0 means use 2-sided differences.
%         dpa(j)=0 (& dpr(j)=0) means hold p(j) fixed.
%         dpa(j)<0 means use 1-sided differences. 
%         If dpa(j)=0, and dpr(j)~=0, relative increments are always used.
% r     = resid(p,P1,...), an array of residuals,
%         initialized by the user before each entry to jnumfit.
% ponoff = column of parameter activities. If ponoff(j)=0, p is held fixed.
% varargin = optional parameters to be passed to resid. If you have
%            passed the arguments A1, A2, etc. to RSLEASQR, then your
%            RESID function should be coded to accept them, i.e.:
%               function [ss,r,f]=resid(p,A1,A2,...).
% ----------- OUTPUT ------------------------------------------------------
% prt = matrix prt(i,j) = dr(i)/dp(j), the Jacobian for use by rsleasqr.
%       where r = the residual vector.  If the original r's are arrays,
%       they will be reshaped to vectors internally for this calculation.
% ------------ GLOBAL FLAG ------------------------------------------------
% The global variable "rsjnumfit_" (note the final underline) signals to
% the world that this function is active, in case there are special
% values that the user wants to preserve during calculation of
% numerical partial derivatives of the fitting parameters.   Upon entry
% to rsjnumfit, rsjnumfit_ is set to 1, and reset to 0 upon return.
% If one wishes to make any use of this variable, the master session 
% should declare it global and set it to 0.
  % RSJNUMFIT                                  page 2/2 
  global rsjnumfit_;                 % Global flag for rsjnumfit.
  rsjnumfit_=1;                      % Signal rsjnumfit in use.
  np=length(p);                      % Dimensions of the problem.
  r=r(:);   nr=size(r,1);            % Reshape r to column.
  ps=p;  prt=zeros(nr,np);           % Save p.  Allocate prt.
  for j=1:np;   p=ps;                % Parameter loop.
    if ~ponoff; continue; end;       % Skip jth loop if inactive.
    meth=sign(dpr(j));               % Method from dpr(j).
    if ~meth,                        % If not from dpr(j),
      meth=sign(dpa(j));             %  method from dpa(j).
      if ~meth;                      % If both = 0,
        continue;  end;  end;        %  skip jth loop.
    meth = 1.5 + 0.5*meth;           % meth = 1 or 2 sided diffs.          
    dp=max(abs(dpr(j)*p(j)),...      % Use relative or
            abs(dpa(j)));            %  absolute increment for p(j).
    if dp==0;                        % If increment=0,
      disp([...                      %  give warning,
        '   Warning: in rsjnumfit, increment ',...
        int2str(j),...               %
	' is zero.']);
	disp('   Parameter, rel. & abs. increments:');
        disp([ps(j),dpr(j),dpa(j)]); % Show parameter & increments.
      continue; end;                 %  & skip jth loop.
    p(j)=ps(j)+dp;                   % Positive increment +dp.
    [ss,r1,f]=feval(...              %
      resid,p,varargin{:});          % Evaluate residuals.
    if ss==-2;                       % Hard abort in resid.
      p(j)=ps(j)-dp;	             % Try negative increment.
      [ss,r1,f]=feval(...            %  
        resid,p,varargin{:});        % Evaluate residuals.
      if ss==-2; continue; end;      % Both sides aborted,
      prt(:,j)=(r-r1(:))/dp;         %  or 1-sided differences.
    elseif meth==1;                  % If 1-sided diff. option,
      prt(:,j)=(r1(:)-r)/dp;         %  use +dp.
    else p(j)=ps(j)-dp;              % Negative increment -dp
      [ss,r2,f]=feval(...            %  for 2-sided diffs.
        resid,p,varargin{:});        % Evaluate residuals.
      if ss==-2;                     % Hard abort in resid.
        prt(:,j)=(r1(:)-r)/dp;       %  Revert to 1-sided diffs.
      else;                          % Both sides evaluated:
        prt(:,j)=(r1(:)-r2(:))...    %  use 2-sided
                 /(2*dp);            %  differences.
      end;   end;   end;             % End ifs & parameter loop.
  rsjnumfit_=0;                      % Signal rsjnumfit ended.