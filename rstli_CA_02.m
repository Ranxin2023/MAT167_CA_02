%% CA2
%Name: Randy Li
%Student ID:917196816
%Date:February 29th 2024

%define m and n
m=50;
n=12;

x = linspace(0, 1, m)';
%create an m*n matrix A
A = ones(m, n); 
for j = 2:n
    A(:, j) = x.^(j-1);  
end

%yi are given by yi = cos ( 6 xi ) i = 1, 2, . . . , m.
y = cos(6*x);

%calculate xhat_normal
xhat_normal = (A' * A) \ (A' * y);

% Display the result
disp('Coefficients of the normal equations result:');
disp(xhat_normal);

% Perform CGS on A
[Q, R] = CGS(A);

% Solve Qy_hat 
y_hat = Q' * y; 

% Solve Rx_hat = y_hat for x_hat
xhat_cgs = R \ y_hat;

% Display the result
disp('Coefficients of the CGS result:');
disp(xhat_cgs);

% Perform Modified Gram-Schmidt on A
[Q, R] = MGS(A);

% Solve Qy_hat = y for y_hat
y_hat = Q' * y;

% Solve Rx_hat = y_hat for x_hat (back substitution)
xhat_mgs = R \ y_hat;  

% Display the result
disp('Coefficients of the MGS result:');
disp(xhat_mgs);

% Perform QR decomposition on A using MATLAB's built-in qr function
[Q, R] = qr(A);

% Solve Qy_hat = y for y_hat
y_hat = Q' * y;  % Since Q is orthogonal, Q'Q = I

% Solve Rx_hat = y_hat for x_hat (back substitution)
xhat_qr = R \ y_hat; 

% Display QR result
disp('Coefficients of the QR result:');
disp(xhat_qr);

% Perform SVD
[U, S, V] = svd(A, 'econ');

% Solve for x_hat using the pseudo-inverse
xhat_svd = V * (S \ (U' * y));

% Display xhat_svd
disp('Coefficients of the polynomial (xhat_svd):');
disp(xhat_svd);

% This one is straightforward with the backslash operator
xhat_backslash = A \ y;

% Display xhat_backslash
disp('Coefficients of the polynomial (xhat_backslash):');
disp(xhat_backslash);

% Combine the results into one matrix, results
Results = [xhat_normal, xhat_cgs, xhat_mgs, xhat_qr, xhat_svd, xhat_backslash];

OUTPUT_DATA_IN_LaTeX_TABLE_FORMAT(Results, '%.16f', 'nomath');
%% Auxiliary function(s)
function [Q,R]=CGS(A)
% Implementation of the Classical Gram-Schmidt Orthogonalization
%  of column vectors of the input matrix A.
%
    [m,n]=size(A);
    V = A;
    Q = zeros(m,n);
    R = zeros(n,n);
    for k=1:n
      for j=1:(k-1)
        R(j,k) = Q(:,j)'*A(:,k);
        V(:,k) = V(:,k)-R(j,k)*Q(:,j);
      end
      R(k,k) = norm(V(:,k));
      Q(:,k) = V(:,k)/R(k,k);
    end
end

%% CALLING SEQUENCE:  [Q,R] = MGS(A);

function [Q,R]=MGS(A)

% Implementation of the Modified Gram-Schmidt Orthogonalization of the 
% column vectors of the input matrix A.

    [m,n]=size(A);
    
    V = A;
    
    Q = zeros(m,n);
    R = zeros(n,n);
    
    for j=1:n
    
      R(j,j) = norm(V(:,j),2);
      Q(:,j) = V(:,j)/R(j,j);
      
      for k=(j+1):n
        R(j,k) = Q(:,j)'*V(:,k);
        V(:,k) = V(:,k)-R(j,k)*Q(:,j);
      end % for(k)
      
    end % for(j)
end

