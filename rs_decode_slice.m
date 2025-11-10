function [decoded_symbols, stats] = rs_decode_slice(codeword_dp, tables)
    % rs_decode_slice - RS(127,121) single-slice decoder (t=3)
    % Input order: [data(1..121), parity(1..6)] (127x1 column vector)

    stats.num_errors = uint8(0);
    stats.uncorrectable = false;

    % Polynomial-descending order for syndrome evaluation:
    % r_126..r_0 = [data(121..1), parity(6..1)]
    
    % +++ BUG FIX 1: horzcat 오류 수정 +++
    % 121x1 열 벡터와 6x1 열 벡터를 각각 전치(transpose)하여
    % 1x121 및 1x6 행 벡터로 만든 후 수평 결합
    cw_desc = [flipud(codeword_dp(1:121)).', flipud(codeword_dp(122:127)).'];
    % 1) Compute syndromes S1..S6
    syndromes = compute_syndromes(cw_desc, tables);

    if all(syndromes == 0)
        decoded_symbols = codeword_dp(1:121); % 121x1 열 벡터 반환
        return;
    end

    % 2) Berlekamp–Massey → lambda, omega
    % +++ BUG FIX 2: poly_utils 인수 전달 +++
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
    % (BUG FIX 2: poly_utils를 이미 인수로 전달받음)
    poly_utils = gf_poly_utils(tables);
    error_values = forney_algorithm_parallel(error_positions, lambda, omega, tables, poly_utils);

    % 5) Apply corrections back to data-parity ordering
    corrected = apply_corrections_dp(codeword_dp, error_positions, error_values, poly_utils);
    decoded_symbols = corrected(1:121); % 121x1 열 벡터 반환
end

% Map degree-indexed error positions to [data(1..121) parity(1..6)] order and correct
% function corrected = apply_corrections_dp(codeword_dp, err_pos, err_val, poly_utils)
%     % codeword_dp는 127x1 열 벡터
%     n = 127;
%     corrected = codeword_dp;
%     for k = 1:length(err_pos)
%         j = double(err_pos(k));           % degree index (0..126)
% 
%         % 다항식 차수 j를 1-based 인덱스 p_desc (1..127)로 변환
%         % j=126 (r_126) -> p_desc=1
%         % j=0   (r_0)   -> p_desc=127
%         p_desc = n - j;
% 
%         % p_desc를 codeword_dp (1..127) 인덱스로 매핑
%         if p_desc <= 121
%             % Data part: p_desc=1 -> idx=121; p_desc=121 -> idx=1
%             idx = 122 - p_desc;
%         else
%             % Parity part: p_desc=122 -> idx=127; p_desc=127 -> idx=122
%             idx = 127 - (p_desc - 122); % 249 - p_desc 와 동일
%         end
% 
%         if idx > 0 && idx <= 127
%             corrected(idx) = poly_utils.gf_add(corrected(idx), err_val(k));
%         end
%     end
% end

% Map degree-indexed error positions to [data(1..121) parity(1..6)] order and correct
function corrected = apply_corrections_dp(codeword_dp, err_pos, err_val, poly_utils)
    % codeword_dp는 LSB-first 127x1 열 벡터입니다.
    % codeword_dp(1) = m0 (degree 6)
    % codeword_dp(121) = m120 (degree 126)
    % codeword_dp(122) = p0 (degree 0)
    % codeword_dp(127) = p5 (degree 5)
    
    n = 127;
    corrected = codeword_dp;
    for k = 1:length(err_pos)
        j = double(err_pos(k));           % degree index (0..126)
        
        % +++ (BUG FIX) New LSB-first index mapping logic +++
        if j >= 6
            % Data part (degree 6..126)
            % j=6   -> idx=1
            % j=126 -> idx=121
            idx = j - 5;
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