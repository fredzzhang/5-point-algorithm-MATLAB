clear all;
close all;
clc;
% -------------------------------------------------------------------------
% Program Introdution:
% This program is a benchmark simulation to verify the 5-point algorithm.

% Prgram flow
% Fix the first camera at the origin and randomly generate a relative pose
% for the second camera. Both cameras are assumed as calibrated, where the 
% intrinsic matrix is an identify matrix. Then randomly generate 5 points 
% in space and project them onto both cameras to obtain matched points in
% two views. 5-point algorithm is then carried out to compute the essential
% matrix and the result is justified.

% Key variables
% thetad - rotation angle in degrees
% theta - rotation angle in radians
% axis - the rotation axis in space, a 3-vector
% R - rotation matrix
% t - translation vector
% p1 - the projection matrix for the first camera
% p2 - the projection matrix for the second camera
% X - points in space, homogeneous coordinate
% x1 - coordinates of matched features in the first view
% x2 - coordinates of matched features in the second view
% eMatrix - the ground truth of essential matrix
% eMatrix5 - the essential matrix computed using 5-point algorithm

% Author: Frederic Zhang
% Last modified: 18 June 2017
% Version: 2.0
% -------------------------------------------------------------------------


% Intrinsic matrix
K = eye(3);

% Rotation angle
thetad = (rand(1) - 0.5) * 360;     % random rotation angle -180~180
theta = thetad * pi / 180;          % rotation angle in radians

% Rotation axis
axis = rand(1, 3);

% Normalize the vector
axis = axis / sqrt(axis(1) ^ 2 + axis(2) ^ 2 + axis(3) ^ 2);    

% Corss product of rotaion vector
v = [0, -axis(3), axis(2); axis(3), 0, -axis(1); -axis(2), axis(1), 0];  
% Rotation matrix according to Rodrigues' Formula
R = eye(3, 3) + sin(theta) * v + (1 - cos(theta)) * v * v;  

% Translation vector
t = [(rand(2, 1) - 0.5) / 10; 1];
% Normalize the last component
t = t / t(3);

% Compute the essential matrix with R and t
tc = [0, -t(3), t(2); t(3), 0, -t(1); -t(2), t(1), 0];
eMatrix = tc * R;
% Normalize the last element
eMatrix = eMatrix / eMatrix(3, 3);

% Projection matrix
p1 = K * [eye(3, 3), [0; 0; 0]];
p2 = K * [R, t];

% 3D points generation
X = rand(4, 5);

% Projection
x1h = p1 * X;
x2h = p2 * X;

% Transformation into inhomogeneous coordinate
x1h = x1h ./ repmat(x1h(3, :), [3, 1]);
x2h = x2h ./ repmat(x2h(3, :), [3, 1]);

x1 = x1h(1:2, :);
x2 = x2h(1:2, :);

% Compute essential matrix using 7-point algorithm
[eMatrix5, num] = fivePoint(x1, x2, K, K);

% Declare a placeholder for reprojection error
err = zeros(1, num);

fprintf('%d solutions are extracted using 5-point algorithm:\n', num);
disp('--------------------------------------------------------------');

% Normalize the last entry of essential matrix and display the results
for i = 1: num
    eMatrix5{i} = eMatrix5{i} / eMatrix5{i}(3, 3);
    fprintf('The %s solution is:\n\n', getOrder(i));
    disp(eMatrix5{i});
    err(i) = sum(diag(x2h' * eMatrix5{i} * x1h));
    fprintf('The associated reprojection error is %f.\n', err(i));
    disp('--------------------------------------------------------------');
end

% Show the ground truth
fprintf('\nThe ground truth of essential matrix is:\n\n');
disp(eMatrix);

% Perform cheirality check to determine the true essential matrix
[E_true, ind] = cheirality(eMatrix5, x1h', x2h', K, K);
fprintf('\nThe true essential matrix is the %s solution:\n\n', getOrder(ind));
disp(E_true);