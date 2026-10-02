% Homework 2, Problem 3
% linearize the 3-state systems at the origin and check the eigenvalues

clear; clc; syms x1 x2 x3 real; % defining the variables

x = [x1; x2; x3]; % state vector
tol = 1e-9;

f1 = [-x1 + x1^2; -x2 + x2^3; x3 - x1^2];

f2 = [-x1; -x1 - x2 - x3 - x1*x3; (x1+1)*x2];

f3 = [-2*x1 + x1^3; -x2 + x1^2; -x3];

F = {f1,f2,f3};

for i = 1:length(F)
    f = F{i};

    % make sure the origin is an equilibrium first
    if any(double(subs(f, x, [0; 0; 0])) ~= 0)
        fprintf('System %d: origin is NOT an equilibrium\n', i);
        continue
    end

    J = jacobian(f, x);
    A = double(subs(J, x, [0; 0; 0]));
    ev = eig(A);
    re = real(ev);

    fprintf('System %d\n', i);
    disp(A);
    fprintf('eigenvalues: %s\n', mat2str(ev.', 4));

    if all(re < -tol)
        fprintf('  -> asymptotically stable at the origin\n\n');
    elseif any(re > tol)
        fprintf('  -> unstable at the origin\n\n');
    else
        fprintf('  -> linearization is inconclusive\n\n');
    end
end
