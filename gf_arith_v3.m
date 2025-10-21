function [tables, g_poly_coeffs] = gf_arithmetic()
% GF(2^7) 연산 테이블 및 생성 다항식 계수 생성
% 
% 출력:
%   tables: GF 테이블 구조체 (exp, log, m)
%   g_poly_coeffs: 생성 다항식 계수 배열

    % --- 표준 파라미터 정의 ---
    m = 7;                          % 필드 차수: GF(2^m)
    p_x_dec = hex2dec('89');        % 원시 다항식 p(x) = x^7 + x^3 + 1 = 0x8B = 139
    
    % --- 1. Generate GF tables ---
    fprintf('1. Start generating GF(2^7) tables...\n');
    tables = generate_gf_tables(m, p_x_dec);
    fprintf('   - gf_exp and gf_log tables generated successfully.\n\n');

    % --- 2. Calculate the coefficients of the generator polynomial g(x) ---
    fprintf('2. Calculate the coefficients of the generator polynomial g(x)...\n');
    g_poly_coeffs = calculate_generator_poly(tables, m);
    
    % --- Result verification ---
    spec_coeffs_hex = {'6D', '22', '64', '44', '40', '7E'};
    fprintf('   - OpenGMSL specification g(x) coefficients (g0~g6): %s\n', strjoin(spec_coeffs_hex, ', '));
    
    calculated_coeffs_hex = arrayfun(@(x) dec2hex(x, 2), g_poly_coeffs, 'UniformOutput', false);
    fprintf('   - 계산된 g(x) 계수 (g0~g6)        : %s\n', strjoin(calculated_coeffs_hex, ', '));

    if isequal(spec_coeffs_hex, calculated_coeffs_hex)
        fprintf('   - 검증 결과: 명세서와 일치합니다! ✅\n\n');
    else
        fprintf('   - 검증 결과: 명세서와 일치하지 않습니다! ❌\n\n');
        fprintf('   - 차이점 분석:\n');
        for i = 1:length(spec_coeffs_hex)
            if ~strcmp(spec_coeffs_hex{i}, calculated_coeffs_hex{i})
                fprintf('     g%d: spec=%s, calc=%s\n', i-1, spec_coeffs_hex{i}, calculated_coeffs_hex{i});
            end
        end
    end

    % --- 3. 시각화 (옵션) ---
    plot_gf_tables = false;  % 디버깅 시에만 true로 설정
    if plot_gf_tables
        fprintf('3. 필드 특성 시각화를 시작합니다...\n');
        plot_exp_table(tables, m);
        plot_multiplication_table(tables, m);
        fprintf('   - 2개의 플롯 생성 완료.\n');
    end
end

% --- 서브 함수: 테이블 생성 ---
function tables = generate_gf_tables(m, p_x_dec)
    % GF(2^m) 지수/로그 테이블 생성
    % 하드웨어: ROM으로 구현
    
    field_size = 2^m;
    gf_exp = zeros(1, field_size - 1);
    gf_log = zeros(1, field_size);
    
    % α^0 = 1
    gf_exp(1) = 1;
    gf_log(1 + 1) = 0; % log(1) = 0
    
    current_val = 1;
    for i = 1:(field_size - 2)
        % α를 곱함 (α=2 이므로 왼쪽 시프트)
        current_val = bitshift(current_val, 1);
        
        % overflow 확인 및 p(x)로 나머지 연산
        if current_val >= field_size
            current_val = bitxor(current_val, p_x_dec);
        end
        
        gf_exp(i + 1) = current_val;
        gf_log(current_val + 1) = i;
    end
    
    tables.exp = gf_exp;
    tables.log = gf_log;
    tables.m = m;
end

% --- 서브 함수: 생성 다항식 계수 계산 ---
function g_poly = calculate_generator_poly(tables, m)
    % g(x) = (x-α^1)(x-α^2)...(x-α^6) 를 전개
    % 하드웨어: 다항식 곱셈기 (루프 언롤링 가능)
    
    % 로컬 헬퍼 함수 정의
    gf_add = @(a, b) bitxor(a, b);

    function result = gf_multiply(a, b)
        if a == 0 || b == 0
            result = 0;
            return;
        end
        log_a = tables.log(a + 1);
        log_b = tables.log(b + 1);
        log_sum = mod(log_a + log_b, 2^m - 1);
        result = tables.exp(log_sum + 1);
    end

    % g(x) 초기값: (x - α^1) = x + α^1
    g_poly = [1, tables.exp(1 + 1)];
    
    for i = 2:6
        root = tables.exp(i + 1); % α^i
        term = [1, root]; % (x - α^i) = x + α^i
        
        % 다항식 곱셈 수행
        new_poly = zeros(1, length(g_poly) + 1);
        for j = 1:length(g_poly)
            for k = 1:length(term)
                prod = gf_multiply(g_poly(j), term(k));
                idx = j + k - 1;
                new_poly(idx) = gf_add(new_poly(idx), prod);
            end
        end
        g_poly = new_poly;
    end
    % 계수 순서를 g0, g1, ... 순으로 맞추기 위해 뒤집음
    g_poly = fliplr(g_poly);
end

% --- 시각화 함수들 (동일) ---
function plot_exp_table(tables, m)
    figure('Name', 'GF(2^7) exponent table visualization');
    
    exponents = 0:(2^m - 2);
    values = tables.exp;
    
    plot(exponents, values, '.-');
    title(['GF(2^', num2str(m), ') exponent table (α^i value)']);
    xlabel('exponent i');
    ylabel('value (10-based)');
    grid on;
    xlim([0, 2^m - 2]);
    ylim([0, 2^m - 1]);
    legend('α^i');
end

function plot_multiplication_table(tables, m)
    figure('Name', 'GF(2^7) multiplication table (16x16)');
    
    function result = gf_multiply(a, b)
        if a == 0 || b == 0
            result = 0;
            return;
        end
        log_a = tables.log(a + 1);
        log_b = tables.log(b + 1);
        log_sum = mod(log_a + log_b, 2^m - 1);
        result = tables.exp(log_sum + 1);
    end

    table_size = 16;
    mult_table = zeros(table_size, table_size);
    
    for i = 0:(table_size - 1)
        for j = 0:(table_size - 1)
            mult_table(i + 1, j + 1) = gf_multiply(i, j);
        end
    end
    
    imagesc(0:(table_size-1), 0:(table_size-1), mult_table);
    colorbar;
    colormap('jet');
    title(['GF(2^', num2str(m), ') multiplication table']);
    xlabel('operand A');
    ylabel('operand B');
    axis square;
    
    % 각 셀에 텍스트 값 추가
    for i = 1:table_size
        for j = 1:table_size
            text(j-1, i-1, num2str(mult_table(i,j)), ...
                 'HorizontalAlignment', 'center', ...
                 'Color', 'w');
        end
    end
end