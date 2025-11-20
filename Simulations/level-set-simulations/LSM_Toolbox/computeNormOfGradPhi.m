function normGradPhi = computeNormOfGradPhi(grid, data, schemeData, signOfV)
%function normGradPhi = computeNormOfGradPhi(grid, data, schemeData, signOfV)

    if(size(signOfV,2) ~= 1)
        error('signOfV must be a scalar in computeNormOfGradPhi');
    end

  normGradPhi = zeros(size(data));
  
  for i = 1 : grid.dim
    
    % Get upwinded derivative approximations.
    [ derivL, derivR ] = feval(schemeData.derivFunc, grid, data, i);
    
    %Figure out the upwind direction
    if(signOfV > 0)
        deriv = derivL;
    else
        deriv = derivR;
    end
    
    % For now, compute the norm squared
    normGradPhi = normGradPhi + deriv.^2;
  end
  
  normGradPhi = sqrt(normGradPhi);
