function error_positions = chien_search(lambda, tables)
    % chien_search - Chien 탐색을 수행하여 오류 위치 다항식의 근을 찾습니다.
    %
    %   이 함수는 RS(127, 121) 코드 (t=3)를 기준으로, σ(x) = 0의 근을
    %   효율적인 하드웨어 구조(병렬 레지스터)로 찾습니다.
    %   σ(α^-j) = σ₀ + σ₁(α⁻ʲ) + σ₂(α⁻ʲ)² + σ₃(α⁻ʲ)³ = 0
    %   이 되는 0-based index 'j' (0..126)를 찾습니다.
    %
    %   입력:
    %       lambda: 1x(t+1) uint8 벡터. BM 알고리즘이 찾은 오류 위치 다항식
    %               σ(x)의 계수 (σ₀, σ₁, σ₂, σ₃)
    %       tables: gf_arithmetic.m에서 생성된 GF 테이블 구조체
    %
    %   출력:
    %       error_positions: 1xV uint8 벡터 (V = 실제 오류 개수).
    %                        오류가 발생한 0-based index (0..126)를 포함.

    % --- 1. 파라미터 및 GF 연산 함수 정의 ---
    m = tables.m;
    n = 2^m - 1; % 127
    max_t = 3;   % RS(127, 121) -> t=3
    num_coeffs = max_t + 1; % σ₀, σ₁, σ₂, σ₃ (총 4개)

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

    % --- 2. Chien 탐색 레지스터 및 상수 초기화 ---
    
    % 하드웨어의 4개 병렬 레지스터 (R₀, R₁, R₂, R₃) 초기화
    % BM이 찾은 lambda 계수를 로드합니다.
    R_regs = zeros(1, num_coeffs, 'uint8');
    len_lambda = length(lambda);
    if len_lambda > num_coeffs
        warning('입력된 다항식 차수가 t=3을 초과합니다.');
        len_lambda = num_coeffs;
    end
    R_regs(1:len_lambda) = lambda(1:len_lambda); % R_regs = [σ₀, σ₁, σ₂, σ₃]
    
    % 하드웨어의 4개 상수 곱셈기 값 (α⁰, α⁻¹, α⁻², α⁻³)
    multipliers = zeros(1, num_coeffs, 'uint8');
    multipliers(1) = 1; % α⁰ = 1
    for i = 1:max_t
        % α⁻ⁱ = exp( (n - i) mod n )
        % (i=1 -> n-1=126), (i=2 -> n-2=125), (i=3 -> n-3=124)
        multiplier_idx = (n - i); % MATLAB 0-based log index
        multipliers(i+1) = tables.exp(multiplier_idx + 1);
    end
    % multipliers = [1, α⁻¹, α⁻², α⁻³]
    
    % 결과를 저장할 빈 배열
    error_positions = [];

    % --- 3. Chien 탐색 127 사이클 반복 ---
    
    % j = 0 부터 126 까지 (총 127번)
    for j = 0:(n-1)
        
        % (a) 계산: Sum = R₀ ⊕ R₁ ⊕ R₂ ⊕ R₃
        %   (j=0일 때) Sum = σ₀ + σ₁ + σ₂ + σ₃ = σ(α⁰)
        %   (j=1일 때) Sum = σ₀ + σ₁(α⁻¹) + σ₂(α⁻²) + σ₃(α⁻³) = σ(α⁻¹)
        %   (j=2일 때) Sum = σ₀ + σ₁(α⁻²) + σ₂(α⁻⁴) + σ₃(α⁻⁶) = σ(α⁻²)
        Sum = uint8(0);
        for i = 1:num_coeffs
            Sum = gf_add(Sum, R_regs(i));
        end
        
        % (b) 판정: Sum이 0이면, j가 오류 위치의 근(역수)입니다.
        if Sum == 0
            % j는 0-based index (0..126)
            error_positions(end + 1) = uint8(j);
        end
        
        % (c) 업데이트: 다음 사이클(j+1)을 위해 레지스터 값 갱신
        % R_i(new) = R_i(old) * (α⁻ⁱ)
        for i = 1:num_coeffs
            R_regs(i) = gf_multiply(R_regs(i), multipliers(i));
        end
    end
end
