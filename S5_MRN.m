% Find equations for yeast glycolysis state variable 1

% define libarary parameters
laurentorder = 0;
polyorder = 2;
usesine = 0;
dyorder = 1;

% search for equation for the 1st state-variable
% We do the ADM method for only one value of lambda 
% see 2nd state varaiable for full sweep.
% clear results for other state variables
clear Theta Thetastring Xi indTheta lambdavec numterms errorv indopt
clear indTheta1 Xi1 numterms1 nT

% pool Data  (i.e., build library of nonlinear time series)
[Theta, Thetastring] = poolDatady(xt,n,polyorder,usesine, laurentorder, dxt(:,5), dyorder);
% %initial lambda value, which is the value used for soft thresholding in ADM

 tol = 7e-10;      %give also good result for upper and lower bound
 lambda = 5e-6;

%  tol = 2e-6;  
%  lambda = 6e-6;
 

jj = 1; % counter
num= 1; % initialize the number of nonzero terms found for the lambda
errorvec= 0;

MaxIter = 1e3;

% for now calculate null space using null function
nT = null(Theta);

[indTheta1, Xi1, numterms1] = ADMinitvary(nT,lambda,MaxIter,tol, plottag);

Thetastring(indTheta1)'
n0Xi = Xi1(Xi1~=0); % terms need to be rearranged to recover coefficients
n0Xi/n0Xi(end)

save('Results/1st_state_variable.mat')