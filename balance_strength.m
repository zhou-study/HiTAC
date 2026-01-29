function B = balance_strength(C)
% 输入:
%   C : K×K 矩阵,  C(i,j) = 簇 i -> 簇 j 的出度
% 输出:
%   B : K×K 对称矩阵, B(i,j) = BalanceStrength(i,j)

    K = size(C,1);
    B = zeros(K);
    eps_val = 0.01;
    % out_deg = sum(C,2);
    % out_deg(out_deg==0) = 1;   % 防止除零
    % C = C ./ out_deg;          % 每行归一化到 [0,1]
    for i = 1:K
        for j = 1:K
            if i ~= j
                out_ij = C(i,j);
                in_ji  = C(j,i);
                num = min(out_ij, in_ji) + eps_val;
                % den = max(out_ij, in_ji) + eps_val;
                % alpha = 0.2;
                B(i,j) = ((2 * num/(out_ij + in_ji + eps_val))) * log(1 + out_ij + in_ji);
                % B(i,j) =  min(out_ij, in_ji) + (out_ij + in_ji) / 2;
                % B(i,j) = out_ij + in_ji + min(out_ij, in_ji);
            end
        end
    end
end

