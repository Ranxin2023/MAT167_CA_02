%% FUNCTION MODIFIED GRAM-SCHMIDT (MGS)
%
% PURPOSE: IMPLEMENTATION OF THE MODIFIED GRAM-SCHMIDT QR FACTORIZATON
% OF THE INPUT MATRIX A.
%
%  Professor Elbridge Gerry Puckett
%  Department of Mathematics
%  University of California, Davis
%  Davis, CA 9516
%
%    EGP Did not write this code. It was presumably written by
%    
%  Professor Naoki Saito
%  Department of Mathematics
%  University of California, Davis
%  Davis, CA 9516
%
%    REVISION HISTORY:
%
%      REVISION 1.00
%
%         Revision 1.00: Sat 19 Feb 2022 01:00:30 PM PST
%         Revision 1.01: 
%
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
