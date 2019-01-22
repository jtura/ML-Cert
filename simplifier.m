%Simplifies a layout removing unnecessary variables
function res = simplifier(layout)
    marked = zeros(1, numel(layout)); %Which subsets are marked for removal
    for i = 1:numel(layout)
        if marked(1,i) == 0 %If it's marked for removal, there's no point in considering it anymore
            for j = 1:numel(layout)
                if marked(1,j)==0
                    if i == j
                    else %We mark for removal all subsets included in the i-th one
                        if size(intersect(layout{i}, layout{j}),2) == size(layout{j},2)
                            marked(1,j) = 1;
                        end
                    end
                end
            end
        end
    end
    res = {};
    for i = 1:numel(layout)
        if marked(1,i)==0
            res = {res{:}, layout{i}};
        end
    end
end