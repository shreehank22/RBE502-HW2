% Homework 2, Problem 5
% cart with inverted pendulum, state z = [x; theta; xdot; thetadot]

clear; clc;
syms m M l g I F th thd xdd thdd x xd real
z = [x; th; xd; thd];

% ---- equations of motion solved for the accelerations ----
eq1 = (m+M)*xdd + m*l*thdd*cos(th) - m*l*thd^2*sin(th) == F;
eq2 = (I + m*l^2)*thdd + m*l*xdd*cos(th) - m*g*l*sin(th) == 0;
acc = solve([eq1, eq2], [xdd, thdd]);
f = [xd; thd; acc.xdd; acc.thdd]; % zdot = f(z, F), nonlinear form (part 2)

% ---- part 1: small angle approximation (cos = 1, sin = theta, thd^2 = 0) ----
D = I*(M+m) + M*m*l^2;
A_sa = [0 0 1 0; 0 0 0 1; 0 -m^2*g*l^2/D 0 0; 0 m*g*l*(M+m)/D 0 0];
B_sa = [0; 0; (I+m*l^2)/D; -m*l/D];

% ---- part 2: Taylor expansion (Jacobians) around theta = 0 ----
vars = [x th xd thd F];
A0 = simplify(subs(jacobian(f, z), vars, [0 0 0 0 0]));
B0 = simplify(subs(diff(f, F),     vars, [0 0 0 0 0]));
disp('A (Taylor, theta = 0)'); disp(A0);
disp('B (Taylor, theta = 0)'); disp(B0);

same_A = all(all(isAlways(simplify(A0 - A_sa) == 0)));
same_B = all(isAlways(simplify(B0 - B_sa) == 0));
fprintf('Taylor matches small angle: A %d, B %d\n', same_A, same_B);

% ---- bonus: same thing around theta = pi ----
Api = simplify(subs(jacobian(f, z), vars, [0 pi 0 0 0]));
Bpi = simplify(subs(diff(f, F),     vars, [0 pi 0 0 0]));
disp('A (theta = pi)'); disp(Api);
disp('B (theta = pi)'); disp(Bpi);

% ---- part 3: numbers and state feedback ----
mv = 0.16; Mv = 0.48; lv = 0.5; gv = 9.8; Iv = mv*lv^2/12;
An = double(subs(A0, [m M l g I], [mv Mv lv gv Iv]));
Bn = double(subs(B0, [m M l g I], [mv Mv lv gv Iv]));
disp('A, B numeric'); disp(An); disp(Bn);
fprintf('open loop eigenvalues: %s\n', mat2str(eig(An).', 4));
fprintf('rank of ctrb = %d\n', rank(ctrb(An, Bn)));

poles = [-2 -3 -4 -5];
K = place(An, Bn, poles);
fprintf('K = %s\n', mat2str(K, 5));
fprintf('closed loop eigenvalues: %s\n', mat2str(eig(An - Bn*K).', 4));

% ---- bonus: symbolic gain ----
% A and B have the form A = [0 0 1 0; 0 0 0 1; 0 a 0 0; 0 b 0 0], B = [0 0 c d]'
syms a b c d s k1 k2 k3 k4 p0 p1 p2 p3
A_s = [0 0 1 0; 0 0 0 1; 0 a 0 0; 0 b 0 0];
B_s = [0; 0; c; d];
cp = expand(det(s*eye(4) - (A_s - B_s*[k1 k2 k3 k4])));
cdes = s^4 + p3*s^3 + p2*s^2 + p1*s + p0; % desired polynomial
Ks = solve(coeffs(cp - cdes, s), [k1 k2 k3 k4]);
disp('symbolic gains in terms of a, b, c, d and the p_i:');
disp(simplify(Ks.k1)); disp(simplify(Ks.k2)); disp(simplify(Ks.k3)); disp(simplify(Ks.k4));

% plug the numbers into the symbolic result and compare with place()
pc = poly(poles); % [1 p3 p2 p1 p0]
vals = [a b c d p3 p2 p1 p0];
nums = [An(3,2) An(4,2) Bn(3) Bn(4) pc(2) pc(3) pc(4) pc(5)];
K_sym = double(subs([Ks.k1 Ks.k2 Ks.k3 Ks.k4], vals, nums));
fprintf('K from symbolic formula = %s\n', mat2str(K_sym, 5));

% ---- nonlinear simulation with F = -K z ----
figure;
for th0 = [0.1 0.3]
    dyn = @(t, q) cartpend(q, -K*q, mv, Mv, lv, gv, Iv);
    [t, Z] = ode45(dyn, [0 5], [0; th0; 0; 0]);
    subplot(2,1,1); plot(t, Z(:,1)); hold on; grid on; ylabel('x [m]');
    subplot(2,1,2); plot(t, Z(:,2)); hold on; grid on; ylabel('\theta [rad]'); xlabel('t [s]');
end
subplot(2,1,1); legend('\theta_0 = 0.1', '\theta_0 = 0.3');

figdir = fullfile(fileparts(mfilename('fullpath')), '..', 'report', 'figures');
if ~exist(figdir, 'dir'), mkdir(figdir); end
exportgraphics(gcf, fullfile(figdir, 'p5.pdf'));


% full nonlinear cart-pendulum dynamics
function dz = cartpend(q, F, m, M, l, g, I)
    th = q(2); thd = q(4);
    Mm  = [m+M, m*l*cos(th); m*l*cos(th), I + m*l^2];
    rhs = [F + m*l*thd^2*sin(th); m*g*l*sin(th)];
    acc = Mm \ rhs;
    dz = [q(3); q(4); acc];
end
