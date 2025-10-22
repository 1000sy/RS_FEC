%function plot_encoder_results(input_bits, encoded_words, crc_enable)
    % Visualizes the results of the rs_encoder function.
    %
    % Inputs:
    %   - input_bits: The original 19-bit/word bitstream fed to the encoder.
    %   - encoded_words: The 128x1 array of 19-bit words from the encoder.
    %   - crc_enable: Flag indicating if CRC was enabled.
    
    %% Plot 1: Data vs Parity Distinction (Stem Plot)
    figure('Name', 'Encoder Output: Data vs Parity');
    word_indices = 0:127;
    
    % Data part (N=0 to 120)
    stem(word_indices(1:121), bitand(encoded_words(1:121), hex2dec('3FFFF')), ...
         'b.', 'MarkerSize', 10, 'DisplayName', 'Data/CRC Words');
    hold on;
    
    % Parity part (N=121 to 126)
    stem(word_indices(122:127), encoded_words(122:127), ...
         'r*', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Parity Words');
    
    % Parity Expansion part (N=127)
    stem(word_indices(128), encoded_words(128), ...
         'k.', 'MarkerSize', 10, 'DisplayName', 'Parity Expansion (Empty)');
    
    grid on;
    hold off;
    title('Encoder Output: Data/CRC and Parity Words', 'FontSize', 14);
    xlabel('Word Index (N)', 'FontSize', 12);
    ylabel('Word Value (18-bit Decimal)', 'FontSize', 12);
    legend('show', 'Location', 'northwest');
    xlim([-5 135]);
    set(gca, 'FontName', 'Malgun Gothic');
    
    
    %% Plot 2: Bit Pattern Comparison (Image Plot)
    figure('Name', 'Encoder Input vs Output Bit Matrix');
    BIT_WIDTH = 19;
    
    % 1. Prepare Input Bit Matrix
    if crc_enable
        num_data_words = 120;
    else
        num_data_words = 121;
    end
    input_matrix_raw = reshape(input_bits, BIT_WIDTH, num_data_words)';
    input_matrix = zeros(128, BIT_WIDTH);
    input_matrix(1:num_data_words, :) = fliplr(input_matrix_raw);
    
    % 2. Prepare Output Bit Matrix
    output_matrix = de2bi(encoded_words, BIT_WIDTH, 'left-msb');
    
    % Plot side-by-side for comparison
    subplot(1, 2, 1);
    imagesc(input_matrix);
    title('Before Encoding: Input Bit Matrix', 'FontSize', 14);
    xlabel('Bit Position (18:is_k, 17:d17..0:d0)', 'FontSize', 12);
    ylabel('Word Index', 'FontSize', 12);
    colormap('gray');
    set(gca, 'FontName', 'Malgun Gothic');
    
    subplot(1, 2, 2);
    imagesc(output_matrix);
    title('After Encoding: Output Bit Matrix', 'FontSize', 14);
    xlabel('Bit Position (18:is_k, 17:d17..0:d0)', 'FontSize', 12);
    ylabel('Word Index', 'FontSize', 12);
    colormap('gray');
    set(gca, 'FontName', 'Malgun Gothic');
    
 %   end
    