function gf = gf_arithmetic()
    % GMSL RS-FEC(127, 121) 사양에 맞는 GF(2^7) 연산 라이브러리
    %
    % 사용법:
    %   gf = gf_arithmetic(); % 라이브러리 초기화
    %   sum_val = gf.add(a, b);
    %   prod_val = gf.multiply(a, b);
    %
    % 반환값 (gf 구조체):
    %   - tables: exp, log 테이블 포함
    %   - add: @(a, b) ... (GF 덧셈 함수 핸들)
    %   - multiply: @(a, b) ... (GF 곱셈 함수 핸들)
    %   - divide: @(a, b) ... (GF 나눗셈 함수 핸들)
    %   - poly_multiply: @(p1, p2) ... (GF 다항식 곱셈 함수 핸들)
    %   - g_poly_coeffs: 검증된 생성 다항식 계수
    
        % --- 표준 파라미터 정의 ---
        m = 7; % 필드 차수: GF(2^m)
        
        % 원시 다항식 p(x) = x^7 + x^3 + 1
        % 이진수: 10001001 -> 0x89 (137)
        % 하드웨어 LFSR에서는 8번째 비트(x^7)가 넘어가면 
        % 나머지 부분인 x^3 + 1 (0x09)와 XOR 연산을 수행합니다.
        p_x_feedback = hex2dec('89'); 
    
        % --- 1. GF 테이블 생성 ---
        gf.tables = generate_gf_tables(m, p_x_feedback);
    
        % --- 2. 연산 함수 핸들 정의 ---
        % 덧셈 (단순 XOR)
        gf.add = @(a, b) bitxor(a, b);
        
        % 곱셈 (룩업 테이블 사용)
        gf.multiply = @(a, b) gf_multiply_lookup(a, b, gf.tables);
        
        % 나눗셈 (룩업 테이블 사용)
        gf.divide = @(a, b) gf_divide_lookup(a, b, gf.tables);
    
        % 다항식 곱셈
        gf.poly_multiply = @(p1, p2) gf_poly_multiply_func(p1, p2, gf);
    
        % --- 3. 생성 다항식 계수 계산 및 검증 ---
        [g_poly_coeffs, is_valid] = verify_generator_poly(gf);
        if is_valid
            fprintf('GF(2^7) 라이브러리가 성공적으로 생성되었으며, 생성 다항식이 명세서와 일치함을 확인했습니다. ✅\n');
            gf.g_poly_coeffs = g_poly_coeffs;
        else
            error('생성된 생성 다항식이 명세서와 일치하지 않습니다! ❌');
        end
    
    end
    
    % --- 서브 함수: 테이블 생성 ---
    function tables = generate_gf_tables(m, p_x_feedback)
        field_size = 2^m;
        
        gf_exp = zeros(1, field_size - 1, 'uint8');
        gf_log = zeros(1, field_size, 'uint8');
        
        % α^0 = 1
        current_val = 1;
        gf_exp(1) = current_val;
        gf_log(current_val + 1) = 0;
        
        for i = 1:(field_size - 2)
            % α를 곱함 (α=2 이므로 왼쪽 시프트)
            current_val = bitshift(current_val, 1);
            
            % overflow 확인 및 p(x)로 나머지 연산
            if current_val >= field_size
                current_val = bitxor(current_val, field_size + p_x_feedback);
            end
            
            gf_exp(i + 1) = current_val;
            gf_log(current_val + 1) = i;
        end
        
        tables.exp = gf_exp;
        tables.log = gf_log;
        tables.m = m;
    end
    
    % --- 서브 함수: GF 곱셈 (테이블 룩업) ---
    function result = gf_multiply_lookup(a, b, tables)
        if a == 0 || b == 0
            result = 0;
            return;
        end
        log_a = tables.log(a + 1);
        log_b = tables.log(b + 1);
        log_sum = mod(double(log_a) + double(log_b), 2^tables.m - 1);
        result = tables.exp(log_sum + 1);
    end
    
    % --- 서브 함수: GF 나눗셈 (테이블 룩업) ---
    function result = gf_divide_lookup(a, b, tables)
        if b == 0
            error('Error: Division by zero in GF');
        end
        if a == 0
            result = 0;
            return;
        end
        log_a = tables.log(a + 1);
        log_b = tables.log(b + 1);
        log_diff = mod(double(log_a) - double(log_b), 2^tables.m - 1);
        result = tables.exp(log_diff + 1);
    end
    
    % --- 서브 함수: GF 다항식 곱셈 ---
    function new_poly = gf_poly_multiply_func(p1, p2, gf)
        len1 = length(p1);
        len2 = length(p2);
        new_poly = zeros(1, len1 + len2 - 1);
        
        for i = 1:len1
            for j = 1:len2
                prod = gf.multiply(p1(i), p2(j));
                idx = i + j - 1;
                new_poly(idx) = gf.add(new_poly(idx), prod);
            end
        end
    end
    
    
    % --- 서브 함수: 생성 다항식 검증 ---
    function [g_poly_coeffs, is_valid] = verify_generator_poly(gf)
        % g(x) = (x-α^1)(x-α^2)...(x-α^6) 를 전개합니다.
        % 명세서에 따르면 j0=1, t=3 이므로 2t=6개의 근을 가집니다.
        
        % g(x) 초기값: (x - α^1) = x + α^1
        g_poly = [1, gf.tables.exp(1 + 1)]; % [1, α^1]
        
        for i = 2:6 % 2t
            root = gf.tables.exp(i + 1); % α^i
            term = [1, root];            % (x - α^i) = x + α^i
            
            g_poly = gf.poly_multiply(g_poly, term);
        end
        
        % 계수 순서를 g0, g1, ... 순으로 맞추기 위해 뒤집음
        g_poly_coeffs = fliplr(g_poly(1:end-1));
    
        % 명세서 값과 비교
        spec_coeffs_hex = {'6D', '22', '64', '44', '40', '7E'};
        spec_coeffs_dec = hex2dec(spec_coeffs_hex)';
        
        is_valid = isequal(g_poly_coeffs, spec_coeffs_dec);
    end