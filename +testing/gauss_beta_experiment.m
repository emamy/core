function gauss_beta_experiment

in = general.interpolation.KernelInterpol;
svm = general.regression.ScalarEpsSVR;
svm.eps = 5;
svm.C = 1000;
k = kernels.GaussKernel;

sv = 1:30;
dims = 1;
fval = repmat([0 5 0 7.5 0 6 -5 0 4 -8],1,3);

sv = repmat(sv,dims,1);
minnorm = Inf;
bestgamma = 0;
minnormsvr = Inf;
bestgammasvr = 0;
for gamma = 1:1:15

    k.Gamma = gamma;
    K = k.evaluate(sv,sv);
    in.K = K;
    a = in.interpolate(fval);

    newnorm = sqrt(a' * K * a);
    if newnorm < minnorm
        minnorm = newnorm;
        bestgamma = gamma;
    end
    
    svm.K = K;
    [ai, b, svidx] = svm.regress(fval);
    
    newnormsvr = sqrt(ai' * K(svidx,svidx) * ai);
    if newnormsvr < minnormsvr
        minnormsvr = newnormsvr;
        bestgammasvr = gamma;
    end

end

x = 1:.1:max(sv(1,:));
x = repmat(x,size(sv,1),1);

k.Gamma = bestgamma;
in.K = k.evaluate(sv,sv);
a = in.interpolate(fval);
fx = a' * k.evaluate(sv,x);
figure;
subplot(1,2,1);
plot(sv,fval,'r*',x,fx,'b');
title(sprintf('Best gamma:%f, ||f||= %2.5f',bestgamma,minnorm));
axis tight;

k.Gamma = bestgammasvr;
svm.K = k.evaluate(sv,sv);
[ai, b, svidx] = svm.regress(fval);
fx = ai' * k.evaluate(sv(:,svidx),x) + b;
subplot(1,2,2);
plot(sv,fval,'r*',x,fx,'b',sv(1,svidx),fval(svidx),'blacks');
title(sprintf('Best gamma:%f, ||f||= %2.5f',bestgammasvr,minnormsvr));
axis tight;

end

