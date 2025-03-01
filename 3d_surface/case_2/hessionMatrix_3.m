function hessian_3 = hessionMatrix_3(a)
    normA = norm(a); % Compute the norm of vector a
    t = a / normA;   % Normalize vector a

    Id3 = eye(3);    % 3x3 Identity matrix

    % Compute the gradient
    grad = (Id3 - (t * t')) / normA;

    % Compute matrix_temp
    matrix_temp = Id3 - (t * t');

    % Compute the Hessian matrix
    hessian_3 = - ( ((grad(:,3) * t') + (t * grad(3,:))) + t(3) * grad ) * normA;
    hessian_3 = hessian_3 + t(3) * matrix_temp + ( t * matrix_temp(3,:) + matrix_temp(:,3) * t' );
    hessian_3 = hessian_3 / (2 * normA * normA);
end