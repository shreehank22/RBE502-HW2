% Homework 2, Problem 1

clear; clc; syms x1 x2 real; % defining the variables

f1 = [-x1 + x1*x2; -x2];

f2 = [-x2-x1*(1 - x1^2 - x2^2);x1 - x2*(1 - x1^2 - x2^2)];

f3 = [ x2*(1-x1^2);-(x1+x2)*(1-x1^2)];

f4 = [-x1-x2;2*x1-x2^3];

F = {f1,f2,f3,f4}; % putting the functions in a cell array

tol = 1e-9;

% state vector
x = [x1; x2];

% verify if the origin is an equilibrium point for each system
for i = 1:length(F)
    f = F{i};
    eq_point = subs(f, x, [0; 0]); % substitute x1=0 and x2=0
    if all(double(eq_point) == 0) % check if the result is a zero vector
        fprintf('System %d has the origin as an equilibrium point.\n', i);
    else
        fprintf('System %d does NOT have the origin as an equilibrium point.\n', i);
    end
end


% evaluate the Jacobian at the origin for each system and store the results
% calculating the eigenvalues of the Jacobian matrices at the origin to determine stability


A = cell(1, length(F)); % preallocate a cell array to store the Jacobians at the origin
eigen_values = cell(1, length(F)); % preallocate a cell array to store the eigenvalues

for i=1:length(F)
    f= F{i};
    J = jacobian(f, x);
    J_at_origin = subs(J, x, [0; 0]);
    A{i} = double(J_at_origin);
    eigenvalues = eig(A{i});
    eigen_values{i} = eigenvalues;
    fprintf('A for System %d:\n', i);
    disp(A{i});
    fprintf('Eigenvalues for System %d:\n', i);
    disp(eigenvalues);
end


% Determine the stability of each system based on the eigenvalues
for i = 1:length(eigen_values)
    real_values = real(eigen_values{i});
    if all(real_values < -tol)
        fprintf('System %d is asymptotically stable at the origin.\n', i);
    elseif any(real_values > tol)
        fprintf('System %d is unstable at the origin.\n', i);
    else
        fprintf('System %d: linearization is inconclusive (some Re(lambda) = 0).\n', i);
    end
end
