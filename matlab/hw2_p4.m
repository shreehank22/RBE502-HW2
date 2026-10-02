% Homework 2, Problem 4
% pole placement for (a) the linear system and (b) the DC motor

clear; clc; syms k1 k2 s;
tol = 1e-9;

% DC motor parameters
R = 2; L = 0.5; Km = 0.1; J = 0.01; b = 0.02;

% system data, A12 of the motor is +K/L, exactly as written in the assignment
A_all = {[1 3; 3 1], [-R/L, Km/L; Km/J, -b/J]};
B_all = {[1; 0], [1/L; 0]};
p_all = {[-1+2j, -1-2j], [-5+5j, -5-5j]}; % desired closed-loop poles
names = {'(a) linear system', '(b) DC motor'};

for i = 1:2
    A = A_all{i}; B = B_all{i}; p = p_all{i};
    fprintf('===== %s =====\n', names{i});

    % 1) stability of the open loop
    ev = eig(A);
    fprintf('eig(A) = %s\n', mat2str(ev.', 4));
    if all(real(ev) < -tol)
        disp('open loop: asymptotically stable');
    elseif any(real(ev) > tol)
        disp('open loop: unstable');
    else
        disp('open loop: marginally stable');
    end

    % 2) controllability
    Co = ctrb(A, B);
    fprintf('rank(ctrb) = %d, det(ctrb) = %g\n', rank(Co), det(Co));

    % 3) gain by matching coefficients of the characteristic polynomial
    Acl = A - B*[k1 k2];
    char_cl = expand(det(s*eye(2) - Acl));
    pd = real(poly(p)); % [1 a1 a0] of the desired polynomial
    char_des = pd(1)*s^2 + pd(2)*s + pd(3);
    eqs = coeffs(char_cl - char_des, s); % each coefficient has to vanish
    sol = solve(eqs, [k1 k2]);
    K_hand = double([sol.k1, sol.k2]);

    % same thing with place() as a check
    K = place(A, B, p);

    fprintf('K by coefficient matching = %s\n', mat2str(K_hand, 6));
    fprintf('K from place()            = %s\n', mat2str(K, 6));
    fprintf('closed-loop poles = %s\n\n', mat2str(eig(A - B*K).', 4));
end
