function [dJdx,dgdx] = gradients(X,p)
%GRADIENTS calculates the gradient of the objective and constraint 
%   Uses sparse finite difference to avoid unnecessary calculations

%     if isMATLABReleaseOlderThan('R2023b')
%         error('This sparse finite difference method requires R2023b or later')
%     end
    % making a dummy user-supplied gradient of zero and comparing
    %fungrad = @(X) deal(simulation(X,p), zeros(numObj,numDV));
    %[valid,err] = checkGradients(@fungrad,X,Display="on")

    % instead of comparing against zero, just call the fd directly
    % this is copied from checkGradients line 136-138
    %[~,JacCineqTrans_fd,JacCeqTrans_fd] = finitedifferences(xPerturb,[],fun, ...
    %    -infBound,infBound,[],cIneq(:),cEq(:),1:sizes.nVar,options,sizes,[],JacCineqTrans_fd, ...
    %    JacCeqTrans_fd,finDiffFlags,[]);

    % syntax above is confusing, use a wrapper that's less confusing
    % https://www.mathworks.com/matlabcentral/answers/1891147-how-does-matlab-compute-numerical-jacobians
  finDiffOpts.DiffMinChange = 0;
  finDiffOpts.DiffMaxChange = inf;
  finDiffOpts.TypicalX = ones(length(X),1);
  finDiffOpts.FinDiffType = 'forward';
  finDiffOpts.GradObj = 'off';
  finDiffOpts.GradConstr = 'off';
  finDiffOpts.UseParallel = false;
  finDiffOpts.FinDiffRelStep = sqrt(eps);
  finDiffOpts.Hessian = [];

  sizes.nVar = length(X);
  sizes.mNonlinIneq = 0;
  sizes.mNonlinEq = 0;
  sizes.xShape = size(X);

  finDiffFlags.fwdFinDiff = true;
  finDiffFlags.scaleObjConstr = false;
  finDiffFlags.chkComplexObj = false;
  finDiffFlags.isGrad = false;
  finDiffFlags.hasLBs = false(sizes.nVar,1);
  finDiffFlags.hasUBs = false(sizes.nVar,1);
  finDiffFlags.chkFunEval = false;
  finDiffFlags.sparseFinDiff = false;

  lenVarIn = 0;
  funValCheck = false;
  gradflag = false;
  hessflag = false;
  funfcn = optimfcnchk(@(x) J_n_wrapper(x, p, n),'fmincon',lenVarIn,funValCheck,gradflag,hessflag);
  J_n = feval(funfcn{3},X);
  dJ_n_dx = sfdnls(X,J_n,[],[],[],funfcn,[],[],finDiffOpts,sizes,finDiffFlags);

    g = feval(funfcn{},X);
    dg_dx = sfdnls(X,g,)

end

function J_n = J_n_wrapper(X,p,n)
    J = simulation(X,p);
    J_n = J(n);
end

