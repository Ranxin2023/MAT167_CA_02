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
