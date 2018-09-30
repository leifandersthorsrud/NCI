function zo = bmffvar_z2day(z, index)

[n, m , t]  = size(z);
T           = numel(index);

zo = zeros(n, m, T);

cnt = t; 
for i = T : - 1 : 1
    zo(:, :, i) = z(:, :, cnt);        
    if index(i) == 1
        cnt = cnt - 1;
    end
    
    if cnt == 0
        %break
        try
            idx = (i - 1 : - 1 : 1);
            zo(:, :, idx) = repmat(zo(:, :, i), [1 1 numel(idx)]);
            break
        catch            
            % do nothing, just break
            break
        end
    end
end