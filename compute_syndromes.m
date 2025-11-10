function syndromes = compute_syndromes(received_symbols, tables)
    % compute_syndromes - RS(127, 121) 코드의 신드롬 6개를 계산합니다.
    %
    %   이 함수는 호너의 방법(Horner's method)을 사용하여 127개의 수신 심볼
    %   다항식 R(x)를 GF(2^7)의 원소 α^1부터 α^6까지에서 평가합니다.
    %   S_i = R(α^i) = ((...(r_126*α^i + r_125)*α^i + ...)*α^i + r_0)
    %
    %   입력:
    %       received_symbols: 1x127 uint8 벡터 (수신된 슬라이스, r_126 ... r_0 순서)
    %       tables:           gf_arithmetic.m에서 생성된 GF 테이블 구조체
    %
    %   출력:
    %       syndromes:        1x6 uint8 벡터 (S1, S2, S3, S4, S5, S6)
    
    % --- 1. 파라미터 및 GF 연산 함수 정의 ---
    
    m = tables.m;
    n = 2^m - 1 ; % 127
    
    % RS(127, 121)은 2t = 6개의 신드롬이 필요합니다.
    num_syndromes = 6; 
    
    if length(received_symbols) ~= n
        error('입력 심볼의 길이가 %d여야 합니다 (현재 %d).', n, length(received_symbols));
    end
    
    % 로컬 GF 연산 함수 (gf_poly_utils.m의 핵심 부분)
    function result = gf_add(a, b)
        result = bitxor(uint8(a), uint8(b));
    end
    
    function result = gf_multiply(a, b)
        if a == 0 || b == 0
            result = 0;
            return;
        end
        log_a = tables.log(a + 1);
        log_b = tables.log(b + 1);
        log_sum = mod(double(log_a) + double(log_b), n);
        result = tables.exp(log_sum + 1);
    end

    % --- 2. 신드롬 계산 준비 ---
    
    % 하드웨어의 6개 병렬 레지스터(누산기)에 해당, 6개의 신드롬 값을 저장하는 레지스터를 준비한다는 말임임
    S_regs = zeros(1, num_syndromes, 'uint8');
    
    % 하드웨어의 6개 상수 곱셈기에 해당 (α^1, α^2, ..., α^6)
    alpha_powers = zeros(1, num_syndromes, 'uint8');
    for i = 1:num_syndromes
        % tables.exp(1) = α^0 = 1
        % tables.exp(2) = α^1
        % ...
        % tables.exp(i + 1) = α^i
        alpha_powers(i) = tables.exp(i + 1); % alpha 의 거듭제곱을 tables.exp 테이블에서 찾아서 저장
    end

    % --- 3. 호너의 방법(Horner's Method) 수행 ---
    %
    % 127 사이클 동안 6개의 신드롬을 병렬로 계산합니다.
    % S_new = (S_old * α^i) + r_k
    %
    % received_symbols(1) = r_126 (x^126의 계수)
    % ...
    % received_symbols(127) = r_0 (x^0의 계수)
    
    for j = 1:n  % 127 사이클 (심볼) 반복 (j=1은 r_126, j=127은 r_0)
        r_k = received_symbols(j);
        
        for i = 1:num_syndromes % 6개 신드롬 병렬 계산
            alpha_i = alpha_powers(i);
            S_old = S_regs(i);
            
            % S_new = (S_old * α^i) + r_k
            S_new = gf_add(gf_multiply(S_old, alpha_i), r_k);
            
            S_regs(i) = S_new;
        end
    end

    syndromes = S_regs;
end
