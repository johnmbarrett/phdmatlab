function h = discreteEntropy(x)
    h = -sum(plogp(epmf(x)));
end