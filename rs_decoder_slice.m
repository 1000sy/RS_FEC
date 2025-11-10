function [decoded_symbols, stats] = rs_decoder_slice(codeword_dp, tables)
    % rs_decode_slice - RS(127,121) single-slice decoder (t=3)
    % Input order: [data(1..121), parity(1..6)]

    stats.num_errors = uint8(0);
    stats.uncorrectable = false;

    % Polynomial-descending order for syndrome evaluation:
    % r_126..r_0 = [data(121..1), parity(6..1)]
    cw_desc = [codeword_dp(121:-1:1), codeword_dp(127:-1:122)];

    % 1) Compute syndromes S1..S6
    syndromes = compute_syndromes(cw_desc, tables);
    if all(syndromes == 0)
        decoded_symbols = codeword_dp(1:121);
        return;
    end

    % 2) Berlekamp–Massey → lambda, omega
    [lambda, omega] = berlekamp_massey(syndromes, tables);

    % 3) Chien search → error degree indices (0..126)
    error_positions = chien_search(lambda, tables);
    if length(error_positions) > 3
        stats.uncorrectable = true;
        decoded_symbols = codeword_dp(1:121);
        return;
    end

    % 4) Forney algorithm (parallel implementation)
    poly_utils = gf_poly_utils(tables);
    error_values = forney_algorithm_parallel(error_positions, lambda, omega, tables, poly_utils);

    % 5) Apply corrections back to data-parity ordering
    corrected = apply_corrections_dp(codeword_dp, error_positions, error_values);
    decoded_symbols = corrected(1:121);
    stats.num_errors = uint8(length(error_positions));
end

% Map degree-indexed error positions to [data(1..121) parity(1..6)] order and correct
function corrected = apply_corrections_dp(codeword_dp, err_pos, err_val)
    n = 127;
    corrected = codeword_dp;
    for k = 1:length(err_pos)
        j = double(err_pos(k));           % degree index (0..126)
        p_desc = n - j;                   % 1-based index in r_126..r_0 vector
        if p_desc <= 121
            idx = 122 - p_desc;           % data index 1..121
        else
            idx = 249 - p_desc;           % parity index 122..127
        end
        corrected(idx) = bitxor(corrected(idx), err_val(k));
    end
end

