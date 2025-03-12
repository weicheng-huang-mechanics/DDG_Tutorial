function hessian_2 = hessionMatrix_2(a)
    normA = norm(a); % Compute the norm of vector a
    t = a / normA;   % Normalize vector a

    Id3 = eye(3);    % 3x3 Identity matrix

    % Compute the gradient
    grad = (Id3 - (t * t')) / normA;

    % Compute matrix_temp
    matrix_temp = Id3 - (t * t');

    % Compute the Hessian matrix
    hessian_2 = - ( ((grad(:,2) * t') + (t * grad(2,:))) + t(2) * grad ) * normA;
    hessian_2 = hessian_2 + t(2) * matrix_temp + ( t * matrix_temp(2,:) + matrix_temp(:,2) * t' );
    hessian_2 = hessian_2 / (2 * normA * normA);
end