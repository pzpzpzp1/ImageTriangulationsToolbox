function approx = Approximator(degree)
    approx = {};
    if degree == 0
        approx.render = @constantRender;
        approx.computeEnergy = @constantComputeEnergy;
        approx.computeGradient = @constantComputeGradient;
    else
        error('Unknown approximator');
    end
end