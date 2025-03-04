function hessian_1 = hessionMatrix_1(a)
    normA = norm(a); % Compute the norm of vector a
    t = a / normA;   % Normalize vector a

    Id3 = eye(3);    % 3x3 Identity matrix

    % Compute the gradient
    grad = (Id3 - (t * t')) / normA;

    % Compute matrix_temp
    matrix_temp = Id3 - (t * t');

    % Compute the Hessian matrix
    hessian_1 = - ( ((grad(:,1) * t') + (t * grad(1,:))) + t(1) * grad ) * normA;
    hessian_1 = hessian_1 + t(1) * matrix_temp + ( t * matrix_temp(1,:) + matrix_temp(:,1) * t' );
    hessian_1 = hessian_1 / (2 * normA * normA);
end