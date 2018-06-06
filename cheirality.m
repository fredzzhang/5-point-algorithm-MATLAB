function [E_true, iE, nmax] = cheirality(E, p, q, K1, K2)
% Function: Return the true essential matrix in a bunch of estimates
% 
% The cheirality check is an operation to eliminate false solutions of
% essential matrices computed using 5-point algorithm or 7-point algorithm
% by reprojecting 3D points back to the associated camera positions. The
% solution with most points in front of both cameras is considered as the
% true solution
% 
% Usage:
%       [E_true, iE, nmax] = cheirality(E)
% where
%       E_true - essential matrix with most points in front of the camera
%       iE - the index of the true essential matrix
%       nmax - the number of points in front of the camera
%       E - several estimates of the essential matrix in a cell array
%       p - point correspondences in the first view
%       q - point correspondences in the second view
%       K1 - the intrinsic matrix of the first camera
%       K2 - the intrinsic matrix of the second camera
% 
% Author: Zhen Zhang
% Institue: Australian National University
% Last modified: 6 Jun. 2018

nE = numel(E);      % get the number of solutions
nmax = 0;            % intialize
iE = 1;

W = [0, -1, 0;...
     1, 0, 0;...
     0, 0, 1];
 
Z = [0, 1, 0;...
     -1, 0, 0;...
     0, 0, 0];

if size(p, 2) ~= 3
    error('Invalid format of points. Reformat to Nx3.');
end
 
if nE == 1
    % If there is only one solution, return it
    E_true = E{1};      
else

for i = 1: nE           % go through all E matrices
    [U, ~, V] = svd(E{i});
    for j = 1: 2        % go through all R and t
        t = U * Z * U';
        t = [t(3, 2); t(1, 3); t(2, 1)];
        switch(j)
            case 1
                R = U * W * V';
            case 2
                R = U * W' * V';
        end
        
        P1 = K1 * [eye(3), zeros(3, 1)];    % camera matrix
        P2 = K2 * [R, t];
        npoint = size(p, 1);
        X = zeros(3, npoint);               % placeholder for 3d points
        for k = 1: npoint
            % get 3d points via triangulation
            X(:, k) = tri3d(p(k, :), q(k, :), P1, P2);
        end
        % project those points back to both cameras
        p1p = P1 * [X; ones(1, npoint)];
        p2p = P2 * [X; ones(1, npoint)];
        if sum((p1p(3, :) > 0) .* (p2p(3, :) > 0)) > nmax
            % note down maximum number of points in front of both cameras
            nmax = sum((p1p(3, :) > 0) .* (p2p(3, :) > 0)); 
            % note down the index
            iE = i;        
        end
    end
end

% Return the true solution
E_true = E{iE}; 

end

end