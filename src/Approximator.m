function approx = Approximator(degree)
    approx = {};
    if degree == 0
        approx.render = @constantRender;
        approx.computeEnergy = @constantComputeEnergy;
        % approx.computeGradient = @constantComputeGradient;
    elseif degree == 1
        approx.computeEnergy = @linearComputeEnergy;
        approx.render = @linearRender;
        % approx.computeGradient = @linearComputeGradient;
    else
        error('Unknown approximator');
    end
end