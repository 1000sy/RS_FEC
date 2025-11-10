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
    data_words_19bit = bi2de(fliplr(reshaped_bits), 'left-msb');
    
    %% 3. Prepare RS Encoder Input Block (Data + CRC)
    rs_input_words_19bit = zeros(121, 1, 'uint32');
    
    if crc_enable
        % --- CRC Enabled Mode ---
        data_for_crc = data_words_19bit(1:120);
        
        % Calculate 18-bit CRC value from 120x19-bit words
        crc_result_18bit = compute_crc_lfsr(data_for_crc);
        
        rs_input_words_19bit(1:120) = data_for_crc;
        % The 121st word is {is_k=0, din=crc_result} as per spec Table 5-7
        rs_input_words_19bit(121) = crc_result_18bit;
        
        stats.crc_value = dec2hex(crc_result_18bit, 5);
    else
        % --- CRC Disabled Mode ---
        rs_input_words_19bit = data_words_19bit; % Use 121 words directly
    end
    
    %% 5. Data Slicing
    % Hardware equivalent: Bit wiring/selection
    % Extract 18-bit data portion {din[17:0]} for RS encoding
    data_part_18bit = bitand(rs_input_words_19bit, hex2dec('3FFFF'));
    
    symbols_slice_A = uint8(bitand(bitshift(data_part_18bit, -12), 63)); % din[17:12]
    symbols_slice_B = uint8(bitand(bitshift(data_part_18bit, -6),  63)); % din[11:6]
    symbols_slice_C = uint8(bitand(data_part_18bit, 63));               % din[5:0]
    
    % Extract is_k vector for final repacking
    is_k_vec = uint8(bitget(rs_input_words_19bit(1:121), BIT_WIDTH));
    
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
    
    %% --- Helper Function: CRC-18 Calculation ---
    function crc_val = compute_crc_lfsr(data_words_19bit)
        % Implements the 18-bit CRC calculation based on Figure 5-14
        % with corrected parameters as per the specification text.
        
        % (B) Initial value: all-ones as per spec
        lfsr = hex2dec('3FFFF'); % 18-bit register initialized to all 1s
        
        % (C) Correct poly_mask for g(x) = x^18 + x^15 + ... + 1
        % using Galois LFSR (left shift, MSB feedback) convention.
        poly_mask = hex2dec('BEA7');
        
        for i = 1:length(data_words_19bit)
            word = data_words_19bit(i);
            
            % (A) Process all 19 bits, LSB-first
            for j = 0:18 % din[0]...din[17], is_k
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
        crc_val = bitand(lfsr, hex2dec('3FFFF')); % Final 18-bit result
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
    
        num_parity = length(g_poly_coeffs); % Should be 6
        lfsr_state = zeros(1, num_parity, 'uint8');
        
        for i = 1:length(data_symbols)
            input_symbol = data_symbols(i);
            feedback = gf_add(input_symbol, lfsr_state(end));
            
            for j = num_parity:-1:2
                term = gf_multiply(feedback, g_poly_coeffs(j));
                lfsr_state(j) = gf_add(lfsr_state(j-1), term);
            end
            lfsr_state(1) = gf_multiply(feedback, g_poly_coeffs(1));
        end
        
        parity_symbols = lfsr_state';
        encoded_slice = [data_symbols; parity_symbols]; % 121 data + 6 parity
    end
    
    % %% --- Helper Function: Repacking Slices ---
    % function final_words = repack_slices(slice_A, slice_B, slice_C, is_k_vec)
    %     % Repacks the three 127-symbol encoded slices into 128 final words.
    %     final_words = zeros(128, 1, 'uint32');
    % 
    %     % 1. Repack Data/CRC part (N = 0 to 120)
    %     data_A = slice_A(1:121);
    %     data_B = slice_B(1:121);
    %     data_C = slice_C(1:121);
    % 
    %     data_part_18bit = bitor(bitor(bitshift(uint32(data_A), 12), ...
    %                                   bitshift(uint32(data_B), 6)), ...
    %                                   uint32(data_C));
    % 
    %     % Combine the 18-bit data with the original 1-bit is_k
    %     final_words(1:121) = bitor(data_part_18bit, bitshift(uint32(is_k_vec), 18));
    % 
    %     % 2. Repack Parity part (N = 121 to 126)
    %     parity_A = slice_A(122:127);
    %     parity_B = slice_B(122:127);
    %     parity_C = slice_C(122:127);
    % 
    %     for i = 1:6
    %         pA = parity_A(i);
    %         pB = parity_B(i);
    %         pC = parity_C(i);
    % 
    %         parity_word_18bit = bitor(bitor(bitshift(uint32(pA), 12), ...
    %                                         bitshift(uint32(pB), 6)), ...
    %                                         uint32(pC));
    % 
    %         % is_k for parity words is always 0
    %         final_words(121 + i) = parity_word_18bit;
    %     end
    % 
    %     % 3. Final word (N=127) is for Parity Expansion, set to 0
    %     final_words(128) = 0;
    % end

%% --- Helper Function: Repacking Slices (CORRECTED) ---
function final_words = repack_slices(slice_A, slice_B, slice_C, is_k_vec)
    % Repacks the three 127-symbol encoded slices into 128 final words.
    final_words = zeros(128, 1, 'uint32');
    
    % 1. Repack Data/CRC part (N = 0 to 120)
    data_A_7bit = slice_A(1:121);
    data_B_7bit = slice_B(1:121);
    data_C_7bit = slice_C(1:121);
    
    % *** 수정된 부분 시작 ***
    % 7-bit 심볼에서 MSB(is_k)를 제거하고 6-bit 데이터만 추출 (Mask = 63 = 0x3F)
    data_A_6bit = bitand(uint32(data_A_7bit), 63);
    data_B_6bit = bitand(uint32(data_B_7bit), 63);
    data_C_6bit = bitand(uint32(data_C_7bit), 63);
    
    % 3개의 6-bit 데이터를 18-bit로 재조립
    data_part_18bit = bitor(bitor(bitshift(data_A_6bit, 12), ...
                                  bitshift(data_B_6bit, 6)), ...
                                  data_C_6bit);
    % *** 수정된 부분 끝 ***
        
    % Combine the 18-bit data with the original 1-bit is_k
    final_words(1:121) = bitor(data_part_18bit, bitshift(uint32(is_k_vec), 18));
    
    % 2. Repack Parity part (N = 121 to 126)
    parity_A = slice_A(122:127); % 패리티 심볼은 is_k=0 이므로 7비트라도 63 이하임
    parity_B = slice_B(122:127);
    parity_C = slice_C(122:127);
    
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
    end
    
    % 3. Final word (N=127) is for Parity Expansion, set to 0
    final_words(128) = 0;
end