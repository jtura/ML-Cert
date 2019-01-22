%Generates a Heisenberg model of n particles and makes the corresponding
%SdP

%layout = [
%0 0 0 0 0 0 0 0 0 0 0;
%1 0 0 0 0 0 0 0 0 0 0;
%0 0 0 0 0 0 0 0 1 0 0;
%0 1 0 0 0 0 0 0 0 0 0;
%0 0 0 1 1 1 0 0 0 0 0]
%generate_example(11,[],[],[],[],[],[],layout)


function res = generate_example(n, Jxx, Jyy, Jzz, Jx, Jy, Jz, layout)
sx = [0 1; 1 0]; sy = sqrt(-1)*[0 -1; 1 0]; sz = [1 0; 0 -1]; s0 = eye(2);
%Generate the terms of the Hamiltonian
H = {};
%1-body terms;
for i = 1:n
    h = Jx(1,i)*sx + Jy(1,i)*sy + Jz(1,i)*sz;
    H{1,i} = h;
end
%2-body terms
for i = 1:n
    h = Jxx(1,i)*kron(sx,sx) + Jyy(1,i)*kron(sy,sy) + Jzz(1,i)*kron(sz,sz);
    H{2,i} = h;
end
res = H;

%We check now if the given layout is compatible with H; i.e., if it is
%sufficient to give a lower bound on the energy
cl = complete_layout(layout);
flag = 1;
cl
for i = 1:2
    for j = 1:n
        if cl(i,j) == 0
            flag = 0;
        end
    end
end
if (flag == 0)
    warning('Layout is too small. Cannot compute objective function\n');
end

%We now assign a pair of indices to every small h term (i.e., the state
%from which to compute the exp value).



return;


%Now we generate a layout at random and solve the corresponding SdP
cvx_begin
    variable Z(n,n) hermitian toeplitz
    dual variable Q
    minimize( norm( Z - P, 'fro' ) )
    Z == hermitian_semidefinite( n ) : Q;
cvx_end
end

function res = complete_layout(layout)
res = layout;
n=size(layout,2);
indices = find(layout'); %find where the ones are
for i = 1:size(indices,1)
    x = 1+floor((indices(i,1)-1)/n);
    y = 1 + mod(indices(i,1)-1,n);
    for fila = 1:x
        for col = y + (0:(x-fila))
            res(fila, mod(col-1,n) + 1) = 1;
        end
    end
end
end