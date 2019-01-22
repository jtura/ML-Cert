%Estimates the size of the SDP, returns the number of free variables.
%Example: ojimetro({[0 1 2], [1 2 3], [3 0]}) is an sdp that
%you can write with 120 variables (and 3 blocks)
function res = ojimetro(layout)
%We implement |A U B| = |A| + |B| - |A\cap B| recursively. We use a
%dynamic programming procedure to speed up things
    c = canonicalize(layout);
    nvars = go(c);
    s = layout; %simplifier(layout);
    aggregated_nblocks = 0;
    for i = 1:numel(s)
        aggregated_nblocks = aggregated_nblocks + 2^(2*size(s{i},2));
    end
    res = [nvars, numel(s), aggregated_nblocks];
end

function res = go(layout)
    persistent M
    if isempty(M)
        M = containers.Map('KeyType','char','ValueType','double');
    end
    slayout = stringify(layout);
    TF = isKey(M, slayout);
    if (TF == false)
        if size(layout,1) == 1
            res = 2^(2*sum(layout))-1;
            M(slayout) = res;
        else
            A = layout(1,:);
            B = layout(2:size(layout,1),:);
            C = zeros(size(B,1),size(B,2));
            for i = 1:size(C,1)
                for j = 1:size(C,2)
                    C(i,j)=B(i,j)*A(1,j);
                end
            end
            A = sortrows(A')';
            B = sortrows(B')';
            C = sortrows(C')';
            res = go(A) + go(B) - go(C);
            M(slayout) = res;
        end
    else
        res = M(slayout);
    end
end

function res = stringify(layout)
    res = '';
    for i = 1:size(layout,1)
        for j = 1:size(layout,2)
            if layout(i,j)==0
                res = strcat(res,'0');
            else
                res = strcat(res,'1');
            end
        end
        res = strcat(res,'n');
    end
end

function res = canonicalize(layout)
    res = zeros(numel(layout),maxel(layout));
    for i = 1:numel(layout)
        for j = 1:size(layout{i},2)
            res(i,1+layout{i}(1,j)) = 1;
        end
    end
    res = sortrows(res')';
end

function res = maxel(layout)
    res = -1;
    for i = 1:numel(layout)
        for j = 1:size(layout{i},2)
            if res < layout{i}(1,j)
                res = layout{i}(1,j);
            end
        end
    end
end