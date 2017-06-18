function [eMatrix, num] = fivePoint(p1, p2, K1, K2)
% -------------------------------------------------------------------------
% Function Introdution:
% Given a set of correspondences between two images and the intrisic matrix
% of the calibrated camera for both views, compute the essential matrix
% associated with the epipolar geometry using five points

% Inputs:
% p1 - the coordinates of matched features in the first image, expressed in
%   inhomogeneous coordinate. The size is of 2xN, where N is the number of 
%   matched points
% p2 - same as above
% K1 - the intrisic matrix of the calibrated camera from the first view
% K2 - the intrisic matrix of the calibrated camera from the second view

% Outputs:
% eMatrix - the computed essential matrix
% num - the number of real roots

% Function Flow
% Apply DLT and construct the 5x9 matrix A
% Extract solutions from the null space of matrix A. It will be written as
%   a linear combination of four components due to fact that the dimension
%   of null space is four. Three coefficients x, y and z are unknowns
% Impose rank constraint and trace constraint. Wrap the ten equations
%   obtained in a 10x20 matrix M
% Perform Gauss-Jordan elimination using partial pivot. The higher order
%   terms will be truncated, leaving a tenth univariate polynomial in 
%   variable z.
% Solve the polynomial using companion matrix and extract ten roots
% Eliminate the complex roots

% Author: Frederic Zhang
% Last modified: 18 June 2017
% Version: 2.0
% -------------------------------------------------------------------------

% 5-point algorithm

% Error checking
[c1, n1] = size(p1);
[c2, n2] = size(p2);
if((c1 ~= 2) || (c2 ~= 2))
    error('Points are not formated with correct number of coordinates.');
end
if((n1 < 5) || (n2 < 5))
    error('There are not enough points to carry out the operation.');
end


% -------------------------------------------------------------------------
% Process inputs
p1 = p1(:, 1: 5);
p2 = p2(:, 1: 5);

p1 = transpose(K1 \ [p1; ones(1, 5)]);
p2 = transpose(K2 \ [p2; ones(1, 5)]);


% Craft matrix A
x1 = p1(:, 1);
y1 = p1(:, 2);
x2 = p2(:, 1);
y2 = p2(:, 2);
A = [x2.*x1, x2.*y1, x2, y2.*x1, y2.*y1, y2, x1, y1, ones(5, 1)];
% Perform SVD
[~, ~, V] = svd(A);
eMatrix1 = [V(1, 6), V(2, 6), V(3, 6); V(4, 6), ...
    V(5, 6), V(6, 6); V(7, 6), V(8, 6), V(9, 6)];
eMatrix2 = [V(1, 7), V(2, 7), V(3, 7); V(4, 7), ...
    V(5, 7), V(6, 7); V(7, 7), V(8, 7), V(9, 7)];
eMatrix3 = [V(1, 8), V(2, 8), V(3, 8); V(4, 8), ...
    V(5, 8), V(6, 8); V(7, 8), V(8, 8), V(9, 8)];
eMatrix4 = [V(1, 9), V(2, 9), V(3, 9); V(4, 9), ...
    V(5, 9), V(6, 9); V(7, 9), V(8, 9), V(9, 9)];

% -------------------------------------------------------------------------
% One equation from rank constraint

% coefficient cellarray
c11 = {eMatrix1(1, 1), eMatrix2(1, 1), eMatrix3(1, 1), eMatrix4(1, 1)};
c12 = {eMatrix1(1, 2), eMatrix2(1, 2), eMatrix3(1, 2), eMatrix4(1, 2)};
c13 = {eMatrix1(1, 3), eMatrix2(1, 3), eMatrix3(1, 3), eMatrix4(1, 3)};
c21 = {eMatrix1(2, 1), eMatrix2(2, 1), eMatrix3(2, 1), eMatrix4(2, 1)};
c22 = {eMatrix1(2, 2), eMatrix2(2, 2), eMatrix3(2, 2), eMatrix4(2, 2)};
c23 = {eMatrix1(2, 3), eMatrix2(2, 3), eMatrix3(2, 3), eMatrix4(2, 3)};
c31 = {eMatrix1(3, 1), eMatrix2(3, 1), eMatrix3(3, 1), eMatrix4(3, 1)};
c32 = {eMatrix1(3, 2), eMatrix2(3, 2), eMatrix3(3, 2), eMatrix4(3, 2)};
c33 = {eMatrix1(3, 3), eMatrix2(3, 3), eMatrix3(3, 3), eMatrix4(3, 3)};

% row1 is a vector containing the 20 coefficients for the first equation
row1 = cell2mat(t4t4t4(c11, c22, c33)) + cell2mat(t4t4t4(c12, c23, c31)) ...
    + cell2mat(t4t4t4(c13, c21, c32)) - cell2mat(t4t4t4(c13, c22, c31)) ...
    - cell2mat(t4t4t4(c12, c21, c33)) - cell2mat(t4t4t4(c11, c23, c32));

% -------------------------------------------------------------------------
% Nine equations from trace constraint

% coefficient cellarray
e1 = {eMatrix1, eMatrix2, eMatrix3, eMatrix4};
e2 = {eMatrix1', eMatrix2', eMatrix3', eMatrix4'};
 
matPart = t4t4t4(e1, e2, e1);
tracePart = t4t10(trace10(t4t4(e1, e2)), e1);

row33 = [];

for i = 1: 20
    % row33 is a cell array containing 20 3x3 matrices
    row33 = [row33, matPart{i} - 0.5 * tracePart{i}];
end

mask = [1; 2; 3; 4; 5; 6; 7; 8; 9];
row9 = [];

for i = 1: 20
    row9 = [row9, row33(mask)];
    row33 = row33(:, 4: end);
end

% -------------------------------------------------------------------------
% Obtain matrix M
M = [row1; row9];

% Rearrange matrix M
colOrder = [1, 7, 2, 4, 3, 11, 8, 14, 5, 12,...
    6, 13, 17, 9, 15, 18, 10, 16, 19, 20];
M = M(:, colOrder);

% Obtain the reduced row enchelon form
rref_M = rref(M);

% Subtraction
eq_k = subtr(rref_M(5, 11: 20), rref_M(6, 11: 20));
eq_l = subtr(rref_M(7, 11: 20), rref_M(8, 11: 20));
eq_m = subtr(rref_M(9, 11: 20), rref_M(10, 11: 20));

% Factorization
B11 = eq_k(1: 4);
B12 = eq_k(5: 8);
B13 = eq_k(9: 13);
B21 = eq_l(1: 4);
B22 = eq_l(5: 8);
B23 = eq_l(9: 13);
B31 = eq_m(1: 4);
B32 = eq_m(5: 8);
B33 = eq_m(9: 13);


% Calculate determinant
p1 = t3t4(B23, B12) - t3t4(B13, B22);
p2 = t3t4(B13, B21) - t3t4(B23, B11);
p3 = t3t3(B11, B22) - t3t3(B12, B21);

% Univariate polynomial of order 10 with z as the variable
detB = t3t7(p1, B31) + t3t7(p2, B32) + t4t6(p3, B33);

% Normalize the coefficient of highest order
if(detB(1))
    detB = detB / detB(1);
end

% Extracting roots of polynomial with companion matrix
Sol_z = eig([zeros(9, 1), eye(9); ...
             - flip(detB(2: end))]);


% Keep the real roots
num = 0;
for i = 1: 10
    if isreal(Sol_z(i)) 
        num = num + 1;
        Sol_z(num) = Sol_z(i);
    end
end
Sol_z = Sol_z(1: num);

% Compute x and y
z6 = [Sol_z .^ 6, Sol_z .^ 5, Sol_z .^ 4, Sol_z .^ 3, Sol_z .^ 2, Sol_z, ...
    ones(num, 1)];
z7 = [Sol_z .^ 7, z6];

p1z = z7 * p1';
p2z = z7 * p2';
p3z = z6 * p3';

Sol_x = p1z ./ p3z;
Sol_y = p2z ./ p3z;

eMatrix = cell(1, num);
for i = 1: num
    
    eMatrix{i} = Sol_x(i) * eMatrix1 + ...
        Sol_y(i) * eMatrix2 + Sol_z(i) * eMatrix3 + eMatrix4;
    
end

function coe10Cell = t4t4(poly1, poly2)
% Given two 4-term polynomials, compute the product
% Inputs: two 4-component coefficient cellarrays
% Output: one 10-component coefficient cellarray

    coe10Cell = {poly1{1} * poly2{1}, ...                            % x^2
                poly1{1} * poly2{2} + poly1{2} * poly2{1}, ...      % x*y
                poly1{1} * poly2{3} + poly1{3} * poly2{1}, ...      % x*z
                poly1{2} * poly2{2}, ...                            % y^2
                poly1{2} * poly2{3} + poly1{3} * poly2{2}, ...      % y*z
                poly1{3} * poly2{3}, ...                            % z^2
                poly1{1} * poly2{4} + poly1{4} * poly2{1}, ...      % x
                poly1{2} * poly2{4} + poly1{4} * poly2{2}, ...      % y
                poly1{3} * poly2{4} + poly1{4} * poly2{3}, ...      % z
                poly1{4} * poly2{4}};                               % 1
   
end

function coe20Cell = t4t10(poly1, poly2)
% Given a 10-term polynomial and a 4-term polynomial, compute the product
% Inputs: one 10-component and a 4-component coefficient cellarrays
% Output: one 20-component coefficient cellarray

    coe20Cell = {poly1{1} * poly2{1}, ...                          % x^3
                poly1{2} * poly2{1} + poly1{1} * poly2{2}, ...    % x^2*y
                poly1{3} * poly2{1} + poly1{1} * poly2{3}, ...    % x^2*z
                poly1{2} * poly2{2} + poly1{4} * poly2{1}, ...    % x*y^2
                poly1{2} * poly2{3} + poly1{3} * poly2{2} + ...   % x*y*z
                poly1{5} * poly2{1}, ...
                poly1{3} * poly2{3} + poly1{6} * poly2{1}, ...    % x*z^2
                poly1{4} * poly2{2}, ...                          % y^3
                poly1{5} * poly2{2} + poly1{4} * poly2{3}, ...    % y^2*z
                poly1{5} * poly2{3} + poly1{6} * poly2{2}, ...    % y*z^2
                poly1{6} * poly2{3}, ...                          % z^3
                poly1{1} * poly2{4} + poly1{7} * poly2{1}, ...    % z^2
                poly1{2} * poly2{4} + poly1{7} * poly2{2} + ...   % x*y
                poly1{8} * poly2{1}, ...                            
                poly1{3} * poly2{4} + poly1{7} * poly2{3} + ...   % x*z
                poly1{9} * poly2{1}, ...
                poly1{4} * poly2{4} + poly1{8} * poly2{2}, ...    % y^2
                poly1{5} * poly2{4} + poly1{8} * poly2{3} + ...   % y*z
                poly1{9} * poly2{2}, ...
                poly1{6} * poly2{4} + poly1{9} * poly2{3}, ...    % z^2
                poly1{7} * poly2{4} + poly1{10} * poly2{1}, ...   % x
                poly1{8} * poly2{4} + poly1{10} * poly2{2}, ...   % y
                poly1{9} * poly2{4} + poly1{10} * poly2{3}, ...   % z
                poly1{10} * poly2{4}};                            % 1
            
end

function coe20Cell = t4t4t4(poly1, poly2, poly3)
% Given three 4-term polynomials, compute the product
% Inputs: three 4-component coefficient cellarrays
% Output: one 20-component coefficient cellarray

    coe20Cell = t4t10(t4t4(poly1, poly2), poly3);
    
end

function coe10Cell = trace10(cell)
% Given ten matrices in a cellarray, compute the trace seperately
% Input: ten matrices in a cellarray
% Output : ten numbers in a cellarray

    coe10Cell = {trace(cell{1}), ...
                 trace(cell{2}), ...
                 trace(cell{3}), ...
                 trace(cell{4}), ...
                 trace(cell{5}), ...
                 trace(cell{6}), ...
                 trace(cell{7}), ...
                 trace(cell{8}), ...
                 trace(cell{9}), ...
                 trace(cell{10})};
end

function coe7Vec = t3t3(poly1, poly2)
% Given two 3-term polynomials, compute the product
% Inputs: two 3-term vectors
% Output: one 7-term vecotr

    coe7Vec = [poly1(1) * poly2(1), ...                             % z6
               poly1(1) * poly2(2) + poly1(2) * poly2(1), ...       % z5
               poly1(2) * poly2(2) + poly1(1) * poly2(3) ...        % z4
               + poly1(3) * poly2(1), ...
               poly1(1) * poly2(4) + poly1(4) * poly2(1) ...        % z3
               + poly1(2) * poly2(3) + poly1(3) * poly2(2), ...      
               poly1(2) * poly2(4) + poly1(4) * poly2(2) ...        % z2
               + poly1(3) * poly2(3), ...
               poly1(3) * poly2(4) + poly1(4) * poly2(3), ...       % z1
               poly1(4) * poly2(4)];                                % 1
           
end

function coe8Vec = t3t4(poly1, poly2)
% Given one 3-term polynomial and a 4-term polynomial, compute the product
% Inputs: one 3-term vector and a 4-term vector
% Output: one 8-term vecotr

    coe8Vec = [poly1(1) * poly2(1), ...                             % z7
               poly1(1) * poly2(2) + poly1(2) * poly2(1), ...       % z6
               poly1(1) * poly2(3) + poly1(2) * poly2(2) ...        % z5
               + poly1(3) * poly2(1), ...                           
               poly1(1) * poly2(4) + poly1(2) * poly2(3) ...        % z4
               + poly1(3) * poly2(2) + poly1(4) * poly2(1), ...
               poly1(2) * poly2(4) + poly1(3) * poly2(3) ...        % z3
               + poly1(4) * poly2(2) + poly1(5) * poly2(1), ...
               poly1(3) * poly2(4) + poly1(4) * poly2(3) ...        % z2
               + poly1(5) * poly2(2), ...
               poly1(4) * poly2(4) + poly1(5) * poly2(3), ...       % z1
               poly1(5) * poly2(4)];                                % 1
           
end

function coe11Vec = t4t6(poly1, poly2)
% Given one 4-term polynomial and a 8-term polynomial, compute the product
% Inputs: one 4-term vector and a 8-term vector
% Output: one 11-term vecotr

    coe11Vec = [poly1(1) * poly2(1), ...                            % z10
                poly1(1) * poly2(2) + poly1(2) * poly2(1), ...      % z9
                poly1(3) * poly2(1) + poly1(2) * poly2(2) ...       % z8
                + poly1(1) * poly2(3), ...
                poly1(4) * poly2(1) + poly1(3) * poly2(2) ...       % z7
                + poly1(2) * poly2(3) + poly1(1) * poly2(4), ...
                poly1(5) * poly2(1) + poly1(4) * poly2(2) ...       % z6
                + poly1(3) * poly2(3) + poly1(2) * poly2(4)...
                + poly1(1) * poly2(5), ...
                poly1(6) * poly2(1) + poly1(5) * poly2(2) ...       % z5
                + poly1(4) * poly2(3) + poly1(3) * poly2(4) ...
                + poly1(2) * poly2(5), ...
                poly1(7) * poly2(1) + poly1(6) * poly2(2) ...       % z4
                + poly1(5) * poly2(3) + poly1(4) * poly2(4) ...
                + poly1(3) * poly2(5), ...
                poly1(7) * poly2(2) + poly1(6) * poly2(3) ...       % z3
                + poly1(5) * poly2(4) + poly1(4) * poly2(5), ...
                poly1(7) * poly2(3) + poly1(6) * poly2(4) ...       % z2
                + poly1(5) * poly2(5), ...
                poly1(7) * poly2(4) + poly1(6) * poly2(5), ...      % z1
                poly1(7) * poly2(5)];                               % 1
            
end

function coe11Vec = t3t7(poly1, poly2)
% Given one 5-term polynomial and a 7-term polynomial, compute the product
% Inputs: one 5-term vector and a 7-term vector
% Output: one 11-term vecotr

    coe11Vec = [poly1(1) * poly2(1), ...                            % z10
                poly1(2) * poly2(1) + poly1(1) * poly2(2), ...      % z9
                poly1(3) * poly2(1) + poly1(2) * poly2(2) ...       % z8
                + poly1(1) * poly2(3), ...
                poly1(4) * poly2(1) + poly1(3) * poly2(2) ...       % z7
                + poly1(2) * poly2(3) + poly1(1) * poly2(4), ...
                poly1(5) * poly2(1) + poly1(4) * poly2(2)...        % z6
                + poly1(3) * poly2(3) + poly1(2) * poly2(4), ...
                poly1(6) * poly2(1) + poly1(5) * poly2(2)...        % z5
                + poly1(4) * poly2(3) + poly1(3) * poly2(4), ...
                poly1(7) * poly2(1) + poly1(6) * poly2(2) ...       % z4
                + poly1(5) * poly2(3) + poly1(4) * poly2(4), ...
                poly1(8) * poly2(1) + poly1(7) * poly2(2) ...       % z3
                + poly1(6) * poly2(3) + poly1(5) * poly2(4), ...
                poly1(8) * poly2(2) + poly1(7) * poly2(3) ...       % z2
                + poly1(6) * poly2(4), ...
                poly1(8) * poly2(3) + poly1(7) * poly2(4), ...      % z1
                poly1(8) * poly2(4)];                               % 1
            
end

function coe13Vec = subtr(vec1, vec2)
% Given two coefficient vectors, calculate vec1 - z*vec2
% Inputs: two 10-component vectors
% Output: one 13-component vector

    v1 = [0, vec1(1: 3), 0, vec1(4: 6), 0, vec1(7: 10)];
    v2 = [vec2(1: 3), 0, vec2(4: 6), 0, vec2(7: 10), 0];
    coe13Vec = v1 - v2;
end

end