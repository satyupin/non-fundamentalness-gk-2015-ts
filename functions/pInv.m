function  InvA  = pInv(A)
% Use pseudo-inverse
temp=eye(size(A,2));
InvA=A\temp;

end

