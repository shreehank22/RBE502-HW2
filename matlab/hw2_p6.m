% Homework 2, Problem 6
% two-wheeled robot, state q = [x y theta nu omega]

clear; clc;
syms x y th nu om k1 k2 real
q = [x; y; th; nu; om];

% closed loop with the given control law
u1 = -k1*nu - x*cos(th) - y*sin(th);
u2 = -k2*om - th;
f_cl = [nu*cos(th); nu*sin(th); om; u1; u2];

% ---- part 1/2: Vdot for V = q'q, with the control law plugged in ----
Vdot = simplify(expand(2*q.'*f_cl));
disp('Vdot ='); disp(Vdot);

% ---- the points (0, y, 0, 0, 0) are all equilibria ----
disp('f_cl on {x = theta = nu = omega = 0} (y free):');
disp(simplify(subs(f_cl, [x th nu om], [0 0 0 0])).');

% ---- linearization of the closed loop at the origin (k1 = k2 = 1) ----
A_cl = double(subs(jacobian(f_cl, q), [q.' k1 k2], [zeros(1,5) 1 1]));
disp('A_cl at the origin:'); disp(A_cl);
fprintf('eigenvalues: %s\n', mat2str(eig(A_cl).', 4));

% ---- simulation ----
K1 = 1; K2 = 1;
dyn = @(t, s) [s(4)*cos(s(3)); s(4)*sin(s(3)); s(5); ...
               -K1*s(4) - s(1)*cos(s(3)) - s(2)*sin(s(3)); ...
               -K2*s(5) - s(3)];

inits = {[1; 1; 0.5; 0; 0], [2; -1; 1; 0.5; 0], [-1; 2; -0.8; 0; 0.3]};
opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);

figure;
subplot(1,2,1); hold on; grid on; xlabel('x'); ylabel('y'); title('o start, * end');
subplot(1,2,2); hold on; grid on; xlabel('t'); ylabel('V = q^T q');
for i = 1:length(inits)
    [t, S] = ode45(dyn, [0 200], inits{i}, opts);
    fprintf('start %s -> end %s\n', mat2str(inits{i}.', 3), mat2str(S(end,:), 4));
    subplot(1,2,1);
    plot(S(:,1), S(:,2)); plot(S(1,1), S(1,2), 'o'); plot(S(end,1), S(end,2), 'k*');
    subplot(1,2,2); plot(t, sum(S.^2, 2));
end

figdir = fullfile(fileparts(mfilename('fullpath')), '..', 'report', 'figures');
if ~exist(figdir, 'dir'), mkdir(figdir); end
exportgraphics(gcf, fullfile(figdir, 'p6.pdf'));
