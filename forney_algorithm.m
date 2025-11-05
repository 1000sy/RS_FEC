function error_values = forney_algorithm(error_positions, lambda, omega, tables)
    % forney_algorithm - Forney 알고리즘을 사용하여 오류 값을 계산합니다.
    %
    %   e_j = Ω(Xj⁻¹) / Λ'(Xj⁻¹)
    %
    %   이 함수는 Chien 탐색으로 찾은 각 오류 위치 j에 대해,
    %   오류 평가 다항식(omega)과 오류 위치 다항식의 도함수(lambda_prime)를
    %   평가하고, 그 값들을 GF 나눗셈하여 실제 오류 값(magnitude)을 계산합니다.
    %
    %   입력:
    %       error_positions: 1xV uint8 벡터. 오류가 발생한 0-based index (0..126).
    %       lambda:          1x(t+1) uint8 벡터. σ(x) 계수.
    %       omega:           1x(t) uint8 벡터. ω(x) 계수.
    %       tables:          gf_arithmetic.m에서 생성된 GF 테이블 구조체.
    %
    %   출력:
    %       error_values:    1xV uint8 벡터. 각 오류 위치에 해당하는 오류 값.
    
    % --- 1. 파라미터 및 GF 연산 함수 정의 ---
    m = tables.m;
    n = 2^m - 1; % 127
    
    error_count = length(error_positions);
    error_values = zeros(1, error_count, 'uint8');
    
    if error_count == 0
        return; % 오류가 없으면 계산할 것도 없음
    end

    % --- 로컬 GF 스칼라 연산 함수 ---
    
    gf_add = @(a, b) bitxor(uint8(a), uint8(b));
    
    function result = gf_multiply(a, b)
        if a == 0 || b == 0, result = 0; return; end
        log_a = tables.log(a + 1);
        log_b = tables.log(b + 1);
        log_sum = mod(double(log_a) + double(log_b), n);
        result = tables.exp(log_sum + 1);
    end
    
    function result = gf_inverse(a)
        if a == 0
            error('GF(2^m)에서 0의 역원은 정의되지 않습니다.');
        end
        log_a = tables.log(a + 1);
        result = tables.exp(n - log_a + 1);
    end
    
    function result = gf_divide(a, b)
        if b == 0
            error('GF(2^m)에서 0으로 나눌 수 없습니다.');
        end
        if a == 0
            result = 0;
            return;
        end
        inv_b = gf_inverse(b);
        result = gf_multiply(a, inv_b);
    end

    % --- 로컬 GF 다항식 연산 함수 ---
    
    % Horner's method로 다항식 P(x)를 x에서 평가
    function result = gf_poly_eval(P, x)
        result = uint8(0);
        % P(x) = P(end)*x^(n-1) + ... + P(2)*x + P(1) (MATLAB 역순)
        % P = [p0, p1, p2, ...] -> P(x) = p0 + p1*x + p2*x^2 ...
        for k = length(P):-1:1
            result = gf_add(gf_multiply(result, x), P(k));
        end
    end
    
    % Λ'(x) (도함수)를 x에서 평가
    % Λ(x) = σ₀ + σ₁x + σ₂x² + σ₃x³
    % Λ'(x) = σ₁ + σ₃x² (홀수 차수 항만 남음) [cite: 297, 3729-3733]
    function result = gf_poly_derivative_eval(P, x)
        result = uint8(0);
        x_inv_sq = gf_multiply(x, x);
        x_inv_pow = uint8(1); % x^0, x^2, x^4 ...
        
        % k=1 -> P(2) (σ₁), k=3 -> P(4) (σ₃)
        for k = 1:2:(length(P)-1)
            coeff = P(k + 1);
            term = gf_multiply(coeff, x_inv_pow);
            result = gf_add(result, term);
            
            x_inv_pow = gf_multiply(x_inv_pow, x_inv_sq);
        end
    end

    % --- 3. 각 오류 위치에 대해 Forney 알고리즘 수행 ---
    
    for i = 1:error_count
        % (a) Chien 탐색이 찾은 오류 위치(인덱스) j
        j = error_positions(i);
        
        % (b) 근(Root) Xj⁻¹ = α⁻ʲ 계산
        if j == 0
            log_index = 0; % α⁰ = 1
        else
            log_index = n - j;
        end
        Xj_inv = tables.exp(log_index + 1);
        
        % (c) 분자 Ω(Xj⁻¹) 계산
        numerator = gf_poly_eval(omega, Xj_inv);
        
        % (d) 분모 Λ'(Xj⁻¹) 계산
        denominator = gf_poly_derivative_eval(lambda, Xj_inv);
        
        % (e) 오류 값 계산 (나눗셈)
        if denominator == 0
            % 이 경우가 발생하면 정정 불가능 오류임
            warning('Forney 알고리즘 분모가 0입니다. (위치: %d)', j);
            error_values(i) = 0; % 오류 값 계산 실패
        else
            error_values(i) = gf_divide(numerator, denominator);
        end
    end
end