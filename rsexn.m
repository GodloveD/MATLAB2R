function [ ss, R, F, C, X ] = rsexn( tau, t, Y, base )
% rsexn: Fitting function for sums of exponentials.
% Given exponential rate constants p, fit Y vs. t to sum of exponentials
%    plus optional base level or base line.
% Coefficients (C) of the exponentials and base function (level or line)
%    are computed by separation of linear parameters.
% -------------- Input Parameters ---------------------------------------
% tau  = column of exponential time constants in increasing order.
% t    = column vector of time.
% Y    = matrix of observed curve(s) y(t) in rows to be fitted by:
%        sum_on_j [ c(q,j)*exp(-p(j)*t) ] + optional base level or
%        function, where q is the index of the row of Y.
% base = scalar base line flag,
%      = 0 if no base level or line,
%      = 1 for base level b1,
%      = 2 for base function, currently b1 + b2/sqrt(t+1).
% ss_ & p_ are global.  They must be initialized by the user
%        to initial sum of squares and parameters before the
%        fit begins.   They are also output.   See below.
% ----------------Output Parameters -------------------------------------
% ss = sum of squares.
% R  = matrix of residuals Y-F.
% F  = matrix of computed curve(s) to be compared to Y.
% C  = matrix of coefficients of exps. in the above fit,
%      with coefficients of each curve in columns.
% X  = matrix of exponentials.
% ssbest_ & taubest_ are best sum of squares and parameters
%      so far, updated whenever the sum of squares improves.
%      This is a safeguard.   In case of an interrupted run,
%      these results will not be destroyed.
   %---------------- Begin Program --------------------------------------- %
   global ssbest_ taubest_;                       % Best ss & p so far.    % Is there a way to remove globals?
   m = length(t);     t = reshape(t,1,m);         % t = row vector.        % t = t';  ??
   n = length(tau);   k = 1./tau(:);              % k = column.            % n is not used; do we need dot with scalar?
                                                                           %
   switch base;                                   % Check base.            %
      case 0;                                     % If base = 0,           %
         X = exp(-k*t);                           %    matrix of exps.     %
      case 1;                                     % If base=1,             %
         X = [exp(-k*t);                          %   exps. +              %
              ones(size(t))];                     %   base level.          %
      case 2;                                     % If base = 2,           %
         X = [exp(-k*t);                          %   exps. +              %
              ones(size(t));                      %   base level +         %
              1./sqrt(t+1) ];                     %   base function.       %
      otherwise error('Bad value of base');       % Error for other        %
   end;                                           %   values of base.      %
                                                                           %
   C = Y/X;   F = C*X;   R = Y-F;                 % Outputs.               % Should the fit be evaluated within this function?  Maybe fit should be output and evaluation should happen outside.  Might remove the need for globals.  
   ss = sum(R(:).^2);                             % Sum of squares.        %
   if ss < ssbest_,                               % If ss has improved,    %
      ssbest_ = ss;                               %   save best ss         %
      taubest_ = tau;                             %   and tau.             %
   end;                                                                    %

