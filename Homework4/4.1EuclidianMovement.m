function [T1, R1, T2, R2, U, V] = TR_from_E(E)
    % This function calculates the possible values for T and R 
    % from the essential matrix
    [U S V] = svd(E);
    if det(U) ~= 1
        U = U * diag([1 1 -1]);
        S = diag([1 1 -1]) * S;
    end
    e2 = U(:,3);
    Rz = [0 -1 0; 1 0 0; 0 0 1];
    R1 = U*Rz'*V';
    T1_hat = U*Rz*S*U';
    T1 = [T1_hat(3,2);T1_hat(1,3);T1_hat(2,1)];
    Rz = [0 1 0; -1 0 0; 0 0 1];
    R2 = U*Rz'*V';
    T2_hat = U*Rz*S*U';
    T2 = [T2_hat(3,2);T2_hat(1,3);T2_hat(2,1)];
end