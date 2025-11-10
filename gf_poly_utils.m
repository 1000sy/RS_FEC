function utils = gf_poly_utils(tables)
    % gf_poly_utils - GF(2^7) 다항식 연산을 위한 유틸리티 함수 팩토리
    %
    %   이 함수는 gf_arithmetic.m에서 생성된 'tables' 구조체를 입력받아,
    %   GF(2^7) 상에서 다항식 연산을 수행하는 함수 핸들이 포함된 
    %   'utils' 구조체를 반환합니다.
    %
    %   입력:
    %       tables: gf_arithmetic.m에서 생성된 구조체 (exp, log, m 필드 포함)
    %
    %   출력:
    %       utils: 다항식 연산 함수 핸들을 포함하는 구조체
    %           .gf_add, .gf_multiply, .gf_inverse, .gf_divide (스칼라 연산)
    %           .add, .conv, .eval, .derivative, .scale, .shift (다항식 연산)

    m = tables.m;
    field_size = 2^m;
    n = field_size - 1; % 127

    % --- 1. 스칼라(Scalar) GF(2^m) 연산 함수 정의 ---
    % 이 함수들은 내부에서 사용되거나 utils 구조체를 통해 외부로 노출됩니다.

    function result = gf_add(a, b)
        % GF(2^m) 덧셈 (Bitwise XOR)
        result = bitxor(uint8(a), uint8(b));
    end

    function result = gf_multiply(a, b)
        % GF(2^m) 곱셈 (Log/Exp 테이블 사용)
        % gf_arithmetic.m에서 사용한 방식과 동일
        if a == 0 || b == 0
            result = 0;
            return;
        end
        % MATLAB 인덱스는 1부터 시작하므로 (값 + 1)
        log_a = tables.log(a + 1);
        log_b = tables.log(b + 1);
        log_sum = mod(double(log_a) + double(log_b), n);
        result = tables.exp(log_sum + 1);
    end

    function result = gf_inverse(a)
        % GF(2^m) 역원 계산
        if a == 0
            error('GF(2^m)에서 0의 역원은 정의되지 않습니다.');
        end
        log_a = tables.log(a + 1);
        log_inv = mod(n - log_a, n);
        result = tables.exp(log_inv + 1);
    end

    function result = gf_divide(a, b)
        % GF(2^m) 나눗셈 (a / b = a * (1/b))
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

    % --- 2. 다항식(Polynomial) GF(2^m) 연산 함수 정의 ---
    % 다항식은 계수의 행 벡터로 표현됩니다 (P(1)=p0, P(2)=p1, ...)

    function Result = poly_add(P1, P2)
        % 다항식 덧셈 (GF(2^m) 계수별 XOR)
        len1 = length(P1);
        len2 = length(P2);
        max_len = max(len1, len2);
        
        % 짧은 다항식을 0으로 패딩
        P1_padded = [P1, zeros(1, max_len - len1, 'uint8')];
        P2_padded = [P2, zeros(1, max_len - len2, 'uint8')];
        
        Result = bitxor(P1_padded, P2_padded);
    end

    function Result = poly_conv(P1, P2)
        % 다항식 곱셈 (Convolution)
        % (P1 * P2)의 결과 다항식 차수는 (deg1 + deg2)
        len1 = length(P1);
        len2 = length(P2);
        
        if len1 == 0 || len2 == 0
            Result = [];
            return;
        end
        
        Result = zeros(1, len1 + len2 - 1, 'uint8');
        
        for i = 1:len1
            for j = 1:len2
                prod = gf_multiply(P1(i), P2(j));
                idx = i + j - 1;
                Result(idx) = gf_add(Result(idx), prod);
            end
        end
    end

    function Result = poly_eval(P, x)
        % 다항식 평가 (Polynomial Evaluation at x)
        % P(x) = P(1) + P(2)*x + P(3)*x^2 + ...
        % Chien 탐색 및 Forney 알고리즘에 사용
        Result = 0;
        x_pow = 1; % x^0
        for i = 1:length(P)
            term = gf_multiply(P(i), x_pow);
            Result = gf_add(Result, term);
            x_pow = gf_multiply(x_pow, x);
        end
    end

    function Result = poly_derivative(P)
        % 다항식의 형식적 미분 (Formal Derivative)
        % Forney 알고리즘에 사용
        % P(x) = p0 + p1*x + p2*x^2 + p3*x^3 + ...
        % P'(x) = p1 + p3*x^2 + p5*x^4 + ...
        if isempty(P)
            Result = [];
            return;
        end
        
        % 홀수 차수 계수 (p1, p3, p5, ...)만 추출
        odd_coeffs = P(2:2:end);
        
        % 짝수 차수 위치 (0, 2, 4, ...)에 배치
        Result = zeros(1, length(P), 'uint8');
        Result(1:2:length(odd_coeffs)*2-1) = odd_coeffs;
        
        % 불필요한 뒤쪽의 0 제거 (MATLAB의 polyval과 유사하게)
        last_nz = find(Result, 1, 'last');
        if isempty(last_nz)
            Result = 0; % 모든 계수가 0이거나 입력이 p0뿐인 경우
        else
            Result = Result(1:last_nz);
        end
    end

    function Result = poly_scale(P, s)
        % 다항식 스케일링 (P * s)
        % 다항식 P의 모든 계수에 스칼라 s를 GF 곱셈
        Result = zeros(1, length(P), 'uint8');
        for i = 1:length(P)
            Result(i) = gf_multiply(P(i), s);
        end
    end

    function Result = poly_shift(P, m)
        % 다항식 시프트 (P * x^m)
        % 다항식 P 앞에 0을 m개 추가
        if isempty(P) || (length(P)==1 && P(1)==0)
            Result = zeros(1, m, 'uint8');
            return;
        end
        Result = [zeros(1, m, 'uint8'), P];
    end


    % --- 3. 반환할 utils 구조체 생성 ---
    utils.gf_add = @gf_add;
    utils.gf_multiply = @gf_multiply;
    utils.gf_inverse = @gf_inverse;
    utils.gf_divide = @gf_divide;
    
    utils.add = @poly_add;
    utils.conv = @poly_conv;
    utils.eval = @poly_eval;
    utils.derivative = @poly_derivative;
    utils.scale = @poly_scale;
    utils.shift = @poly_shift;
end
