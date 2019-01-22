%if A = a1a2a3 and B = b1b2b3b4, where the mask is 0110100, with ordering [2 3 0] we get
%b1a2a3b2a1b3b4; i.e., Alice has been reshuffled according to the ordering
%given
function res = indexmerge2(A,B,mask,ordering)
[sorted, perm] = sort(ordering);
n1 = sum(mask);
n2 = size(mask,2)-n1;
a = d2b(A,n1);
%Permute a
aa(perm) = a; %We can puit the index in the lhs!! wow
b = d2b(B,n2);
ca = 1; cb = 1; res = 0;
for i = 1:size(mask,2)
    res = res * 2;
    if mask(1,i)==1
        res = res + aa(1,ca);
        ca = ca+1;
    else
        res = res + b(1,cb);
        cb=cb+1;
    end
end
end

function res = d2b(b,k)
res = zeros(1,k);
for i = 1:k
    res(1,k-i+1) = mod(b,2);
    b=b-mod(b,2);
    b=b/2;
end
end