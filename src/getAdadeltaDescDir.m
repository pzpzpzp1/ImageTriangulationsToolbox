% Custom implementation of adadelta and rmsprop. 
% Differs from standard approach in that it has a learning rate per (x,y) pair. So vertices have independent update rates, but x and y are scaled together.
% a.la. https://ruder.io/optimizing-gradient-descent/index.html#gradientdescentoptimizationalgorithms
function [descDir, rates] = getAdadeltaDescDir(grad, rateflag)
    nX = size(grad,1);
    persistent expG2 expTh2;
    if numel(expG2)==0 || size(grad,1)~=size(expG2,1);
        descDir = -grad;
        expG2 = sum(grad.^2,2);
        expTh2 = sum(descDir.^2,2);
        rates = ones(nX,1);
        return;
    end
    
    gamma = .7; 
    eps = 1e-8;
    
    % update decaying average of previous gradient square norms
    expG2 = gamma*expG2 + (1-gamma)*sum(grad.^2,2);
    
    if rateflag == 1
        % AdaDelta
        rates = sqrt((expTh2+eps)./(expG2+eps));
    elseif rateflag == 2
        % RMSprop
    	rates = sqrt(1./(expG2+eps));
    end
    descDir = -rates.*grad;
    
    % update decaying average of previous descent dir square norms
    expTh2 = gamma*expTh2 + (1-gamma)*sum(descDir.^2,2);
end