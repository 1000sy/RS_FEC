function [decoded_bits, stats] = rs_decoder_main(encoded_words, tables, crc_enable)
    % rs_decoder_main - Block-level RS-FEC decoder orchestrator
    % encoded_words: 128x1 uint32, each word is 19-bit {is_k, din[17:0]}

    if numel(encoded_words) ~= 128
        error('encoded_words must have 128 entries');
    end

    % 1) Unpack three slices and is_k vector
    [slice_A, slice_B, slice_C, is_k_vec] = unpack_slices(encoded_words);

    % 2) Per-slice decode
    [decoded_A, stats_A] = rs_decode_slice(slice_A, tables);
    [decoded_B, stats_B] = rs_decode_slice(slice_B, tables);
    [decoded_C, stats_C] = rs_decode_slice(slice_C, tables);

    % 3) Repack corrected 121 symbols
    words_19bit = repack_data(decoded_A, decoded_B, decoded_C, is_k_vec);

    % 4) CRC check (optional)
    if crc_enable
        stats.crc_pass = crc_check(words_19bit);
    else
        stats.crc_pass = true;
    end

    % 5) Convert words to bits (LSB-first per word)
    decoded_bits = words_to_bits(words_19bit);

    % 6) Aggregate stats
    stats.corrected_errors = double(stats_A.num_errors) + double(stats_B.num_errors) + double(stats_C.num_errors);
    stats.uncorrectable = logical(stats_A.uncorrectable || stats_B.uncorrectable || stats_C.uncorrectable);
end

% ===== Helpers =====
function [sliceA, sliceB, sliceC, is_k] = unpack_slices(words)
    % Reverse of encoder repack logic
    BIT_MASK_18 = uint32(hex2dec('3FFFF'));
    sliceA = zeros(127,1,'uint8');
    sliceB = zeros(127,1,'uint8');
    sliceC = zeros(127,1,'uint8');
    is_k   = zeros(121,1,'uint8');

    % Data/CRC part N=1..121
    for n = 1:121
        w = uint32(words(n));
        is_k(n) = uint8(bitget(w,19));
        d = bitand(w, BIT_MASK_18);
        sliceA(n) = uint8(bitand(bitshift(d,-12), uint32(63)));
        sliceB(n) = uint8(bitand(bitshift(d,-6),  uint32(63)));
        sliceC(n) = uint8(bitand(d,                      uint32(63)));
    end

    % Parity part N=122..127 (6 words)
    for i = 1:6
        w = uint32(words(121+i));
        d = bitand(w, BIT_MASK_18);
        idx = 121 + i; % 122..127
        sliceA(idx) = uint8(bitand(bitshift(d,-12), uint32(63)));
        sliceB(idx) = uint8(bitand(bitshift(d,-6),  uint32(63)));
        sliceC(idx) = uint8(bitand(d,                      uint32(63)));
    end
end

function words = repack_data(A, B, C, is_k)
    % Pack 121 corrected data symbols back to 19-bit words
    words = zeros(121,1,'uint32');
    d18 = bitor(bitor(bitshift(uint32(A(1:121)),12), bitshift(uint32(B(1:121)),6)), uint32(C(1:121)));
    words(1:121) = bitor(d18, bitshift(uint32(is_k),18));
end

function pass = crc_check(words)
    % Compute CRC-18 over first 120 words and compare to word 121 (lower 18 bits)
    if numel(words) < 121, pass=false; return; end
    crc_calc = compute_crc_lfsr(words(1:120));
    crc_word = bitand(uint32(words(121)), uint32(hex2dec('3FFFF')));
    pass = (crc_calc == crc_word);
end

function crc_val = compute_crc_lfsr(data_words_19bit)
    % Same bit-serial LFSR as encoder (LSB-first)
    lfsr = uint32(hex2dec('3FFFF'));           % 18-bit all ones
    poly_mask = uint32(hex2dec('BEA7'));       % g(x)=x^18+x^15+...+1
    for i = 1:length(data_words_19bit)
        word = uint32(data_words_19bit(i));
        for j = 0:18
            inb = bitget(word, j+1);
            fb  = bitget(lfsr, 18);
            lfsr = bitshift(lfsr,1);
            lfsr = bitset(lfsr,1, bitxor(inb, fb));
            if fb
                lfsr = bitxor(lfsr, poly_mask);
            end
        end
    end
    crc_val = bitand(lfsr, uint32(hex2dec('3FFFF')));
end

function bits = words_to_bits(words)
    % Flatten 19-bit words to a bit vector (LSB-first per word)
    N = numel(words);
    bits = zeros(N*19,1,'uint8');
    idx = 1;
    for i=1:N
        w = uint32(words(i));
        for b = 0:18
            bits(idx) = uint8(bitget(w, b+1));
            idx = idx + 1;
        end
    end
end

