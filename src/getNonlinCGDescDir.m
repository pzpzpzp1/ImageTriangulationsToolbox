% implements 4 different methods of nonlinear conjugate gradient a.la.wikipidia.
% wikipedia recommends P.R.
function [descDir, beta] = getNonlinCGDescDir(grad, betaFlag)
    persistent prevGrad prevDescDir
    if numel(prevGrad)==0
        prevGrad = zeros(size(grad));
        prevDescDir = prevGrad;
    end
    
    if betaFlag == 1
        beta = norm(grad(:))^2/norm(prevGrad(:))^2; % Fletcher-Reeves
    elseif betaFlag == 2
        beta = grad(:)'*(grad(:)-prevGrad(:))/norm(prevGrad(:))^2; % Polak-Ribiere
    elseif betaFlag == 3
        beta = grad(:)'*(grad(:)-prevGrad(:))/(prevDescDir(:)'*(grad(:)-prevGrad(:))); % Hestenes-Stiefel
    elseif betaFlag == 4
        beta = norm(grad(:))^2/(prevDescDir(:)'*(grad(:)-prevGrad(:))); % Dai-Yuan
    else
        error('unhandled beta flag in nonlinear conjugate gradient.');
    end
    
    % no negative momentum. will either go along the previous descent direction or not be influenced. 
    if isinf(beta) || isnan(beta) || beta < 0
        beta = 0; 
    end
    
    descDir = -grad + beta*prevDescDir;
    prevDescDir = descDir;
    prevGrad = grad;
end