function [res, En] = go1D(layout)
sx = [0 1; 1 0]; sy = sqrt(-1)*[0 -1; 1 0]; sz = [1 0; 0 -1]; s0 = eye(2);
n = size(layout,2);
k = size(layout,1);
n
%min(eig(H5(-0.3,1,1)))
%pause;
%We define the Hamiltonian
H = {};
for i = 0:n-1
    H = {H{:}, {[i], -0.3*sz}};
    H = {H{:}, {[i, mod(i+1,n)], kron(sx,sx)+kron(sy,sy)}};
end
F= [];
energy = 0;
%We convert the layout in a set of masks, one for each variable
V={};
for i = 1:k
    for j = 0:n-1
        if layout(i,1+j)==1
            tmp = zeros(1,0);
            for a = 0:i-1
                tmp = [tmp(1,:), mod(j+a, n)]; %This is the ordering of the state!, so it can be e.g. Charlie-Alice
            end
            V = {V{:}, {tmp,sdpvar(2^i, 2^i)}};
            F = [F, V{numel(V)}{2}>=0, trace(V{numel(V)}{2})==1];
        end
    end
end
%We try now to build the cost function
flag = 0;
for i = 1:numel(H)
    found = 0;
    for j = 1:numel(V)
        if found == 0
            if size(intersect(H{i}{1}, V{j}{1}),2)==size(H{i}{1},2)
                found = 1;
                %Example: H acts on parties 501 and rho on parties 45012.
                %Then we need to compute sum_i,j,a1a2a3,b1b2b3 of
                %(H)_{a1a2a3,b1b2b3} * (rho)_{i|b1b2b3|j,i|a1a2a3|j}
                %fprintf('Joining...');
                %H{i}{1}
                %V{j}{1}
                
                %RIGHT NOW WE JUST ASSUME THAT EVERY TERM TO COMPUTE THE
                %ENERGY IS IN V, AND THEY APPEAR IN THE RIGHT ORDER
                % I.E. FOR 2-BODY N-N SHOULD NOT BE A PROBLEM, BUT MAYBE IF
                % YOU TRY SOMETHING MORE EXOTIC
                energy = energy + trace(H{i}{2} * V{j}{2}); %This has to be fixed, order, etc. matters
            end
        end
    end
    if found == 0
        flag = 1;
    end
end
if flag == 1
    error 'Cannot compute energy with current layout';
end

%Now we put compatibility constraints
for i = 1:numel(V)
    for j = (i+1):numel(V)
        common = intersect(V{i}{1}, V{j}{1}); %indices are ordered here.
        if size(common,2)>0
            %fprintf('New case!!! %d %d', i, j);
            %V{i}{1}
            %V{j}{1}
            %common
            maskLeft = ismember(V{i}{1}, common).*1;
            maskRight = ismember(V{j}{1}, common).*1;
            %maskLeft
            %maskRight
            orderL = intersect(V{i}{1},common, 'stable');
            orderR = intersect(V{j}{1},common, 'stable'); %Order matters now!!
            %For each element of the common RDM
            for a = 0:(2^size(common,2)-1)
                for b = 0:(2^size(common,2)-1)
                    %Now we have an equality constraint for the element
                    %(a,b) of the intersection between the two. We have to
                    %complete it and sum the row-columns with (a,b) fixed.
                    %Note that we also have to undo the permutation; i.e.,
                    %if we have rho on parties [2 3 4 0] and we want to
                    %reduce it to rho on parties [0 2 3] we don't want to
                    %end up with rho on parties [2 3 0]
                    constr = 0;
                    for c = 0:(2^(size(V{i}{1},2)-size(common,2))-1) %Completes the index for the first state
                        %[1+indexmerge(a,c,maskLeft),1+indexmerge(b,c,maskLeft)]
                        %[a b c maskLeft]
                        constr = constr + V{i}{2}(1+indexmerge2(a,c,maskLeft, orderL),1+indexmerge2(b,c,maskLeft, orderL));
                    end
                    for d = 0:(2^(size(V{j}{1},2)-size(common,2))-1) %Completes the index for the second state
                        %[1+indexmerge(a,d,maskRight),1+indexmerge(b,d,maskRight)]
                        %[a b d maskRight]
                        constr = constr - V{j}{2}(1+indexmerge2(a,d,maskRight, orderR),1+indexmerge2(b,d,maskRight, orderR));
                    end
                    F=[F,constr==0];
                end
            end
        end
    end
end

%Run the SdP
optimize(F, energy)

% saving energy as double
En = double(energy);

res = V;
end

function lol = H5(z,x,y)
sx = [0 1; 1 0]; sy = sqrt(-1)*[0 -1; 1 0]; sz = [1 0; 0 -1]; s0 = eye(2);
lol = z*(k5(sz,s0,s0,s0,s0)+k5(s0,sz,s0,s0,s0)+k5(s0,s0,sz,s0,s0)+k5(s0,s0,s0,sz,s0)+k5(s0,s0,s0,s0,sz));
lol = lol + x*(k5(sx,sx,s0,s0,s0)+k5(s0,sx,sx,s0,s0)+k5(s0,s0,sx,sx,s0)+k5(s0,s0,s0,sx,sx)+k5(sx,s0,s0,s0,sx));
lol = lol + y*(k5(sy,sy,s0,s0,s0)+k5(s0,sy,sy,s0,s0)+k5(s0,s0,sy,sy,s0)+k5(s0,s0,s0,sy,sy)+k5(sy,s0,s0,s0,sy));
end

function res = k5(a,b,c,d,e)
res = kron(kron4(a,b,c,d),e);
end

%FailCase
% layout = kron([1 1 1 1 1], [1 1 0 0 0]')
% l(5,1)=1
% go1D(l)
% l = layout;
% l(5,1)=1
% go1D(l)
% l = layout;
% l(4,3)=1
% l(4,1)=1
% go1D(l)