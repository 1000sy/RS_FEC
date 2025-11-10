function error_values = forney_algorithm_parallel(error_positions, lambda, omega, tables, poly_utils)
    % forney_algorithm - Forney 알고리즘을 t=3 병렬 하드웨어 구조로 모델링합니다.
    %
    %   e_j = Ω(Xj⁻¹) / Λ'(Xj⁻¹)
    %
    %   이 함수는 Chien 탐색으로 찾은 최대 3개의 오류 위치에 대해
    %   3개의 병렬 계산 경로를 사용하여 오류 값을 동시에 계산합니다.
    %
    %   입력:
    %       error_positions: 1xV uint8 벡터 (V <= 3). 오류가 발생한 0-based index.
    %       lambda:          1x(t+1) uint8 벡터. σ(x) 계수.
    %       omega:           1x(t) uint8 벡터. ω(x) 계수. (BM에서 계산됨)
    %       tables:          gf_arithmetic.m에서 생성된 GF 테이블 구조체.
    %       poly_utils:      gf_poly_utils.m에서 생성된 유틸리티 구조체.
    %
    %   출력:
    %       error_values:    1xV uint8 벡터. 각 오류 위치에 해당하는 오류 값.
    
    m = tables.m;
    n = 2^m - 1; % 127
    
    error_count = length(error_positions);
    
    if error_count == 0
        error_values = [];
        return; % 오류가 없으면 계산할 것도 없음
    end

    % --- 1. 공통 계산 (Common Logic) ---
    % 하드웨어의 3개 병렬 경로가 모두 공유하는 계산입니다.
    % Λ'(x) (도함수) 다항식 계수 계산
    lambda_prime_poly = poly_utils.derivative(lambda);

    % --- 2. 병렬 계산 경로 초기화 ---
    % 하드웨어의 3개 병렬 계산 경로를 모델링합니다.
    % MATLAB에서는 순차적으로 실행되지만, 각 'if' 블록은
    % RTL에서 error_count >= 1'을 Enable 신호로 갖는
    % 독립적인 하드웨어 블록을 의미합니다.
    
    val1 = uint8(0);
    val2 = uint8(0);
    val3 = uint8(0);

    % --- 3. 병렬 경로 1 (오류 1개 이상일 때 활성화) ---
    if error_count >= 1
        j = error_positions(1);
        
        % (b) 근(Root) Xj⁻¹ = α⁻ʲ 계산
        if j == 0, log_index = 0; else, log_index = n - j; end
        Xj_inv = tables.exp(log_index + 1);
        
        % (c) 분자 Ω(Xj⁻¹) 계산
        num1 = poly_utils.eval(omega, Xj_inv);
        
        % (d) 분모 Λ'(Xj⁻¹) 계산
        den1 = poly_utils.eval(lambda_prime_poly, Xj_inv);
        
        % (e) 오류 값 계산 (나눗셈)
        if den1 ~= 0    % 분모가 0이 아니면 오류 값 계산
            val1 = poly_utils.gf_divide(num1, den1);
        else
            warning('Forney 알고리즘 분모(1)가 0입니다. (위치: %d)', j); % 분모가 0이면 정정 불가능 오류
            val1 = uint8(0); % 정정 불가능 오류
        end
    end
    
    % --- 4. 병렬 경로 2 (오류 2개 이상일 때 활성화) ---
    if error_count >= 2
        j = error_positions(2);
        
        if j == 0, log_index = 0; else, log_index = n - j; end
        Xj_inv = tables.exp(log_index + 1);
        
        num2 = poly_utils.eval(omega, Xj_inv);
        den2 = poly_utils.eval(lambda_prime_poly, Xj_inv);
        
        if den2 ~= 0
            val2 = poly_utils.gf_divide(num2, den2);
        else
            warning('Forney 알고리즘 분모(2)가 0입니다. (위치: %d)', j);
            val2 = uint8(0); % 정정 불가능 오류
        end
    end
    
    % --- 5. 병렬 경로 3 (오류 3개일 때 활성화) ---
    if error_count >= 3
        j = error_positions(3);
        
        if j == 0, log_index = 0; else, log_index = n - j; end
        Xj_inv = tables.exp(log_index + 1);
        
        num3 = poly_utils.eval(omega, Xj_inv);
        den3 = poly_utils.eval(lambda_prime_poly, Xj_inv);
        
        if den3 ~= 0
            val3 = poly_utils.gf_divide(num3, den3);
        else
            warning('Forney 알고리즘 분모(3)가 0입니다. (위치: %d)', j);
            val3 = uint8(0); % 정정 불가능 오류
        end
    end
    
    % --- 6. 최종 결과 취합 (Mux) ---
    % 하드웨어에서 3개의 경로에서 나온 값을 MUX로 선택하는 과정입니다.
    all_values = [val1, val2, val3];
    error_values = all_values(1:error_count);
    
end