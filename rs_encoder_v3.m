%function [encoded_words, stats] = rs_encoder_v3(tables, g_poly_coeffs)
    % GMSL RS-FEC(127, 121) Encoder (Hardware-friendly Model)
    %
    % This model processes 19-bit input words ({is_k, din[17:0]}) and
    % generates 128 output words including data, CRC, and RS parity.
    % It strictly adheres to the hardware implementation principles.
    
    %% 1. Initialization and Input Validation
    crc_enable = 0;
    NUM_WORDS = 121;

    BIT_WIDTH = 19; % Input word width including is_k bit
    input_bits = randi([0 1], NUM_WORDS * BIT_WIDTH, 1);
    is_k_indices = 1 : BIT_WIDTH : (NUM_WORDS * BIT_WIDTH); % is_k bit indices
    input_bits(is_k_indices) = 0;
    
    stats.crc_enabled = crc_enable;
    stats.input_bit_count = length(input_bits);
    
    
    if crc_enable
        expected_bits = 120 * BIT_WIDTH;
    else
        expected_bits = 121 * BIT_WIDTH;
    end
    
    if stats.input_bit_count ~= expected_bits
        error('Input bit count does not match the selected mode (requires %d-bit words).', BIT_WIDTH);
    end
    
    %% 2. Bitstream to 19-bit Word Conversion
    % Hardware equivalent: Serial-to-parallel converter
    reshaped_bits = reshape(input_bits, BIT_WIDTH, [])';
    % 'left-msb' 옵션: 첫 번째 열(is_k)이 MSB, 마지막 열(din[0])이 LSB
    data_words_19bit = bi2de(reshaped_bits, 'left-msb');
    
    %% 3. Extract is_k and 18-bit data separately
    % is_k는 별도로 보관 (CRC/RS 인코더 계산에서 분리)
    is_k_vec_raw = uint8(bitget(data_words_19bit(1:121), BIT_WIDTH)); % is_k 비트 추출
    data_words_18bit = bitand(data_words_19bit, hex2dec('3FFFF'));    % 18-bit 데이터만 추출
    
    %% 4. Prepare RS Encoder Input Block (Data + CRC)
    rs_input_words_18bit = zeros(121, 1, 'uint32');  % 18-bit 데이터만 저장
    rs_input_words_19bit = zeros(121, 1, 'uint32');  % 최종 19-bit 워드 (repack용)
    
    if crc_enable
        % --- CRC Enabled Mode ---
        % CRC는 18-bit 데이터만 처리 (is_k 제외)
        data_for_crc = data_words_18bit(1:120);
        
        % Calculate 18-bit CRC value from 120x18-bit data words
        crc_result_18bit = compute_crc_lfsr(data_for_crc);
        
        rs_input_words_18bit(1:120) = data_for_crc;
        % The 121st word is {is_k=0, din=crc_result} as per spec Table 5-7
        rs_input_words_18bit(121) = crc_result_18bit;
        
        % 19-bit 워드 재구성 (is_k는 원본 유지, 121번째만 is_k=0)
        rs_input_words_19bit(1:120) = bitor(data_words_18bit(1:120), bitshift(uint32(is_k_vec_raw(1:120)), 18));
        rs_input_words_19bit(121) = crc_result_18bit; % is_k = 0
        
        stats.crc_value = dec2hex(crc_result_18bit, 5);
    else
        % --- CRC Disabled Mode ---
        rs_input_words_18bit = data_words_18bit; % 18-bit 데이터만 사용
        rs_input_words_19bit = data_words_19bit; % 19-bit 워드 그대로 사용
    end
    
    %% 5. Data Slicing
    % Hardware equivalent: Bit wiring/selection
    % RS 인코더는 18-bit 데이터만 사용 (is_k는 별도 처리)
    data_part_18bit = rs_input_words_18bit; % 이미 18-bit만 추출됨
    % is_k 벡터는 위에서 추출한 is_k_vec_raw 사용 (CRC 모드에서도 원본 유지)
    if crc_enable
        % CRC 모드: 121번째 워드의 is_k는 항상 0
        is_k_vec = [is_k_vec_raw(1:120), uint8(0)]; % 121번째만 0으로 설정
    else
        is_k_vec = is_k_vec_raw;
    end

    % 6-bit LSB 부분 추출
    lsbs_A = uint8(bitand(bitshift(data_part_18bit, -12), 63)); % din[17:12]
    lsbs_B = uint8(bitand(bitshift(data_part_18bit, -6),  63)); % din[11:6]
    lsbs_C = uint8(bitand(data_part_18bit, 63));               % din[5:0]

    % 7-bit 심볼 재조립: MSB (bit 6) = is_k
    msb = bitshift(uint32(is_k_vec), 6);
    
    symbols_slice_A = bitor(uint32(lsbs_A), msb); % [is_k, din[17:12]]
    symbols_slice_B = bitor(uint32(lsbs_B), msb); % [is_k, din[11:6]]
    symbols_slice_C = bitor(uint32(lsbs_C), msb); % [is_k, din[5:0]]
    
    % uint8로 변환
    symbols_slice_A = uint8(symbols_slice_A);
    symbols_slice_B = uint8(symbols_slice_B);
    symbols_slice_C = uint8(symbols_slice_C);
    
    %% 6. Parallel RS Encoding for each slice
    % Hardware equivalent: Three parallel RS encoder LFSRs (Fig 5-13)
    encoded_slice_A = rs_encode_slice(symbols_slice_A, tables, g_poly_coeffs);
    encoded_slice_B = rs_encode_slice(symbols_slice_B, tables, g_poly_coeffs);
    encoded_slice_C = rs_encode_slice(symbols_slice_C, tables, g_poly_coeffs);
    
    %% 7. Repacking
    % Hardware equivalent: Multiplexers and bit concatenation logic
    encoded_words = repack_slices(encoded_slice_A, encoded_slice_B, encoded_slice_C, is_k_vec);
    
    fprintf('RS-FEC encoding completed successfully.\n');
    
    %end
    
%% --- Helper Function: CRC-18 Calculation (BUG FIXED) ---
function crc_val = compute_crc_lfsr(data_words_18bit)
    % CRC-18 계산: 18-bit 데이터만 처리 (is_k 비트 제외)
    % 입력: 18-bit 데이터 워드 배열 (is_k 없음)
    
    lfsr = uint32(hex2dec('3FFFF')); % 18-bit register initialized to all 1s
    poly_mask = uint32(hex2dec('BEA7'));
    
    for i = 1:length(data_words_18bit)
        word = uint32(data_words_18bit(i)); % 18-bit 데이터 워드
        
        % Process 18 bits only (din[0]...din[17]), LSB-first
        for j = 0:17 % din[0]...din[17] (is_k 제외)
            input_bit = bitget(word, j + 1);
            
            % Hardware: Feedback is taken from the MSB (s17) of the register
            feedback_bit = bitget(lfsr, 18);
            
            % Hardware: Register shifts left by 1 bit
            lfsr = bitshift(lfsr, 1);
            
            % Hardware: Input bit is XORed with feedback and injected at LSB (s0)
            lfsr = bitset(lfsr, 1, bitxor(input_bit, feedback_bit));
            
            % Hardware: If feedback is 1, XOR the register with the poly_mask
            if (feedback_bit)
                lfsr = bitxor(lfsr, poly_mask);
            end
        end
    end
    crc_val = bitand(lfsr, uint32(hex2dec('3FFFF'))); % Final 18-bit result
end



%% --- Helper Function: Single Slice RS Encoding ---
function encoded_slice = rs_encode_slice(data_symbols, tables, g_poly_coeffs)
    % Implements the systematic RS encoder using a Galois LFSR (Fig 5-13)
    
    gf_add = @(a, b) bitxor(a, b);
    function result = gf_multiply(a, b)
        if a == 0 || b == 0, result = 0; return; end
        log_a = tables.log(a + 1);
        log_b = tables.log(b + 1);
        log_sum = mod(double(log_a) + double(log_b), 127);
        result = tables.exp(log_sum + 1);
    end

    % --- (BUG FIX) ---
    % num_parity = length(g_poly_coeffs); % This is 7. WRONG.
    num_parity = 6; % RS(127, 121) has 6 parity symbols
    
    g_coeffs = g_poly_coeffs(1:num_parity); % Extract g0..g5 only
    
    lfsr_state = zeros(1, num_parity, 'uint8'); % 6 registers
    
    for i = 1:length(data_symbols)
        input_symbol = data_symbols(i);
        feedback = gf_add(input_symbol, lfsr_state(end)); % lfsr_state(6)
        
        for j = num_parity:-1:2
            % term = gf_multiply(feedback, g_poly_coeffs(j)); % WRONG
            term = gf_multiply(feedback, g_coeffs(j)); % Use g_coeffs(2..6) -> g1..g5
            lfsr_state(j) = gf_add(lfsr_state(j-1), term);
        end

        lfsr_state(1) = gf_multiply(feedback, g_coeffs(1)); % Use g_coeffs(1) -> g0
    end
    
    % 패리티 심볼 추출: LFSR state를 그대로 사용 (fliplr 없음)
    % lfsr_state(1) = p_0, lfsr_state(6) = p_5
    % RS 코드워드: [data(1..121), parity(1..6)] = [data, p_0, p_1, ..., p_5]
    
    % parity_symbols = lfsr_state';
    % encoded_slice = [data_symbols; parity_symbols]; % 121 data + 6 parity
    parity_symbols = fliplr(lfsr_state);   % [p5 p4 ... p0]
    encoded_slice = [data_symbols; parity_symbols'];
end

%% --- Helper Function: Repacking Slices (CORRECTED) ---
function final_words = repack_slices(slice_A, slice_B, slice_C, is_k_vec)
    % Repacks the three 127-symbol encoded slices into 128 final words.
    final_words = zeros(128, 1, 'uint32');
    
    % 1. Repack Data/CRC part (N = 0 to 120)
    data_A_7bit = slice_A(1:121);
    data_B_7bit = slice_B(1:121);
    data_C_7bit = slice_C(1:121);
    
    % 7-bit 심볼에서 MSB(is_k)를 제거하고 6-bit 데이터만 추출 (Mask = 63 = 0x3F)
    data_A_6bit = bitand(uint32(data_A_7bit), 63);
    data_B_6bit = bitand(uint32(data_B_7bit), 63);
    data_C_6bit = bitand(uint32(data_C_7bit), 63);
    
    % 3개의 6-bit 데이터를 18-bit로 재조립
    data_part_18bit = bitor(bitor(bitshift(data_A_6bit, 12), ...
                                  bitshift(data_B_6bit, 6)), ...
                                  data_C_6bit);
        
    % Combine the 18-bit data with the original 1-bit is_k
    final_words(1:121) = bitor(data_part_18bit, bitshift(uint32(is_k_vec), 18));
    
    % 2. Repack Parity part (N = 121 to 126)
    parity_A = slice_A(122:127); % 패리티 심볼은 is_k=0 이므로 7비트라도 63 이하임
    parity_B = slice_B(122:127);
    parity_C = slice_C(122:127);
    
    parity_expansion = uint32(0);

    for i = 1:6
        pA = parity_A(i);
        pB = parity_B(i);
        pC = parity_C(i);
        
        % 패리티는 7비트(최대 127)이지만, MSB(is_k)는 항상 0이므로
        % 6비트 마스킹(bitand(..., 63))은 사실 필요 없으나,
        % 데이터 경로와 일관성을 맞추기 위해 6비트만 사용한다고 가정합니다.
        % (만약 패리티도 7비트 심볼을 그대로 쓴다면 18비트가 아닌 21비트가 필요합니다)
        % RS 인코더 출력이 7비트이므로, 6비트만 잘라서 쓰는 것이 맞습니다.
        parity_word_18bit = bitor(bitor(bitshift(bitand(uint32(pA), 63), 12), ...
                                        bitshift(bitand(uint32(pB), 63), 6)), ...
                                        bitand(uint32(pC), 63));
        
        % is_k for parity words is always 0 (as per spec table)
        final_words(121 + i) = parity_word_18bit;

        % +++ 상위 1비트(MSB)는 Parity Expansion 워드에 모으기 +++
        % Spec Table 5-7의 Expansion 구조에 맞춰 비트 배치
        % (여기서는 단순화를 위해 순서대로 배치한다고 가정, 실제 Spec에 따라 시프트 조정 필요)
        % 예: pA_msb는 0,3,6... pB_msb는 1,4,7... 위치 등 (여기선 단순 패킹)
        msb_A = bitget(pA, 7); 
        msb_B = bitget(pB, 7);
        msb_C = bitget(pC, 7);

        % pA[i]의 MSB를 expansion 워드의 3*(i-1) 위치에, pB는 +1, pC는 +2 위치에 배치
        parity_expansion = bitset(parity_expansion, 3*(i-1) + 1, msb_A);
        parity_expansion = bitset(parity_expansion, 3*(i-1) + 2, msb_B);
        parity_expansion = bitset(parity_expansion, 3*(i-1) + 3, msb_C);

    end
    
    % 3. Final word (N=127) is for Parity Expansion, set to 0
    final_words(128) = parity_expansion;
end