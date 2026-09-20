function [df, testVariance, factorVariance, chi2Lower, chi2Upper] = chi2Test(V, Q, numberOfObservations, numberOfParams)

% X^(2) test for variance
df = numberOfObservations - numberOfParams;
testVariance = (transpose(V) * Q * V);
factorVariance = testVariance/df;

% testing the factor_variance
alpha = 0.05;
chi2Lower = chi2inv(alpha/2, df);
chi2Upper = chi2inv(1 - alpha/2, df);

if testVariance > chi2Lower && testVariance < chi2Upper
    disp('Model passes the chi-squared test.');
else
    disp('Model fails the chi-squared test.');
end

if testVariance < chi2Lower
    disp('Model failed the test on the lower side')
elseif testVariance > chi2Upper
    disp('Model failed the test on the upper side')
end

end
