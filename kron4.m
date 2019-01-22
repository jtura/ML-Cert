function res = kron4(a,b,c,d)
res = kron(kron(a,b),kron(c,d));
end
%kron4(sx,sz,sx,s0)+kron4(s0,sx,sz,sx)+kron4(sx,s0,sx,sz)+kron4(sz,sx,s0,sx)-2*(kron4(sx,sz,sy,s0)+kron4(s0,sx,sz,sy)+kron4(sy,s0,sx,sz)+kron4(sz,sy,s0,sx)) - (kron4(sy,sz,sy,s0)+kron4(s0,sy,sz,sy)+kron4(sy,s0,sy,sz)+kron4(sz,sy,s0,sy)) +eye(16)*8