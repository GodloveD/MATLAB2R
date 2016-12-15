vmax=1;  km=1;  p=[vmax;km];             % Simulation study. First set
x=(0:.1:2)';    y=zeros(21,1);           % p to known values, then
[s,r,f]=hyperb(p,x,y);                   % generate artificial data.
niter=50; stol=.001; erra=1;             % Request for error analysis.
y=f+.01*randn(21,1);                     % Perturb p.  See how well p
dpr=ones(2,1); dpa=dpr;                  % is recovered by fitting to
ponoff = [];
[ p,se,dep,  r,f,s,  iter,cnvrg,qi ]= ...      % the noisy data.  Note: x & y
    rsleasqr('hyperb','dhyperb', ...           % are arguments to be passed
    p,dpr,dpa,ponoff, ...                 % to hyperb & dhyperb.
    stol,niter,erra,x,y);
disp('Parameters Std. errs. Dependencies');
disp([p,se,dep]);                        % Show parameters et al.
figure
plot(x,y,'o')
hold on
plot(x,f,'k');                   % Plot y and fitted f.
