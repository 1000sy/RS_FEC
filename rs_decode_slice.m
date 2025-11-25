function [decoded_symbols, stats] = rs_decode_slice(codeword_dp, tables)
    % rs_decode_slice - RS(127,121) single-slice decoder (t=3)
    % Input order: [data(1..121), parity(1..6)] (127x1 column vector)

    stats.num_errors = uint8(0);
    stats.uncorrectable = false;

    % Polynomial-descending order for syndrome evaluation:
    % r_126..r_0 = [data(121..1), parity(6..1)]

    %cw_desc = [codeword_dp(121:-1:1); codeword_dp(127:-1:122)];
    cw_desc = [codeword_dp(1:121); codeword_dp(122:127)];
    % cw_desc = [fliplr(codeword_dp(1:121)); fliplr(codeword_dp(122:127))];
    % cw_desc = fliplr(cw_desc(:)')';
    % cw_desc = codeword_dp(127:-1:1);

    % 1) Compute syndromes S1..S6
    syndromes = compute_syndromes(cw_desc, tables);

    if all(syndromes == 0)
        decoded_symbols = codeword_dp(1:121); % 121x1 열 벡터 반환
        return;
    end

    % 2) Berlekamp–Massey → lambda, omega
    [lambda, omega] = berlekamp_massey(syndromes, tables);

    % 3) Chien search → error degree indices (0..126)
    error_positions = chien_search(lambda, tables);
    
    stats.num_errors = uint8(length(error_positions));
    lambda_degree = length(lambda) - 1;

    % 정정 가능 여부
    if stats.num_errors > 3 || stats.num_errors ~= lambda_degree || lambda_degree > 3
        stats.uncorrectable = true;
        decoded_symbols = codeword_dp(1:121);
        return;
    end
    
    if stats.num_errors == 0
        % 신드롬은 0이 아니었으나 BM/Chien 결과 오류가 0개인 경우 (이론상 가능)
        decoded_symbols = codeword_dp(1:121);
        return;
    end

    % 4) Forney algorithm (parallel implementation)
    poly_utils = gf_poly_utils(tables);
    error_values = forney_algorithm_parallel(error_positions, lambda, omega, tables, poly_utils);

    % 5) Apply corrections back to data-parity ordering
    corrected = apply_corrections_dp(codeword_dp, error_positions, error_values, poly_utils);
    decoded_symbols = corrected(1:121); % 121x1 열 벡터 반환
end

function corrected = apply_corrections_dp(codeword_dp, err_pos, err_val, poly_utils)
    %n = 127;
    corrected = codeword_dp;
    for k = 1:length(err_pos)
        j = double(err_pos(k));           % degree index (0..126)

        % +++ (BUG FIX) New LSB-first index mapping logic +++
        if j >= 6
            % Data part (degree 6..126)
            % j=6   -> idx=1
            % j=126 -> idx=121
            % idx = j - 5; 이게 틀린거고
            idx = 127 - j;
        else
            % Parity part (degree 0..5)
            % j=0   -> idx=122
            % j=5   -> idx=127
            idx = 122 + j;
        end
        % +++ (END BUG FIX) +++

        if idx > 0 && idx <= 127
            corrected(idx) = poly_utils.gf_add(corrected(idx), err_val(k));
        else
            warning('Error correction index out of bounds: j=%d, idx=%d', j, idx);
        end
    end
end