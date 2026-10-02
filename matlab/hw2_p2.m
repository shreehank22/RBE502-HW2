% Homework 2, Problem 2
% find all the equilibria of each system and check their stability

clear; clc; syms x1 x2 real; % defining the variables

x = [x1; x2]; % state vector
tol = 1e-9;

f1 = [x1-x1^3+x2;3*x1-x2];

f2 = [(x1-2)*(4+x1-2*x2); x2*(x1-x2)];

f3 = [x2*(3-x1-2*x2); x1*(2-x1-x2)];

F = {f1,f2,f3}; % systems 1-3 are polynomial so solve() works fine

for i = 1:length(F)
    f = F{i};
    sol = solve(f == 0, [x1 x2]); % all real solutions of f(x) = 0
    xe = double([sol.x1, sol.x2]); % one equilibrium per row
    J = jacobian(f, x);

    fprintf('System %d has %d equilibria\n', i, size(xe,1));
    for k = 1:size(xe,1)
        A = double(subs(J, x, xe(k,:).')); % Jacobian at this equilibrium
        ev = eig(A);
        fprintf('xe = (%g, %g), A = %s\n', xe(k,1), xe(k,2), mat2str(A,4));
        fprintf('eig(A) = %s -> %s\n', mat2str(ev.',4), classify_eq(ev, tol));
    end
end


% system 4 has tan terms so I do not trust solve() here
% subtracting the two equations gives x1 = x2 = x (worked out in the report),
% so the equilibria are the roots of g(x) = x - 0.5*tan(pi*x/2) on (-1,1)
f4 = [-0.5*tan(pi*x1/2) + x2; x1 - 0.5*tan(pi*x2/2)];
J4 = jacobian(f4, x);

g = @(t) t - 0.5*tan(pi*t/2);
xs = linspace(-0.99, 0.99, 2000); % even number of points so x = 0 is not on the grid
roots_g = [];
for k = 1:length(xs)-1
    if g(xs(k))*g(xs(k+1)) < 0 % sign change means a root in between
        roots_g(end+1) = fzero(g, [xs(k) xs(k+1)]); %#ok<SAGROW>
    end
end

fprintf('System 4 has %d equilibria in |x| < 1\n', length(roots_g));
for k = 1:length(roots_g)
    r = roots_g(k);
    res = norm(double(subs(f4, x, [r; r]))); % should be ~0
    A = double(subs(J4, x, [r; r]));
    ev = eig(A);
    fprintf('xe = (%.4f, %.4f), residual = %.1e, A = %s\n', r, r, res, mat2str(A,4));
    fprintf('eig(A) = %s -> %s\n', mat2str(ev.',4), classify_eq(ev, tol));
end


% classify an equilibrium from the eigenvalues of its Jacobian
function s = classify_eq(ev, tol)
    re = real(ev); im = imag(ev);
    if all(re < -tol)
        if all(abs(im) < tol), s = 'stable node'; else, s = 'stable focus'; end
    elseif all(re > tol)
        if all(abs(im) < tol), s = 'unstable node'; else, s = 'unstable focus'; end
    elseif any(re > tol) && any(re < -tol)
        s = 'saddle (unstable)';
    else
        s = 'inconclusive (Re(lambda) = 0)';
    end
end
