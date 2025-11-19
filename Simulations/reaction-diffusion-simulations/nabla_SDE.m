function nablaX = nabla_SDE(X)
    nablaX = circshift(X,1)+circshift(X,-1)-2*X;
end