function [lambda, omega] = berlekamp_massey(syndromes, tables)
    % berlekamp_massey - 개선된 Berlekamp-Massey 알고리즘 (Inversionless)
    %
    %   학위 논문 5.2.3절의 'Improved Berlekamp-Massey Algorithm' [cite: 1308-1327]을 구현합니다.
    %   이 함수는 2t번의 반복을 통해 GF 나눗셈 없이 Key Equation을 풉니다.
    %
    %   입력:
    %       syndromes: 1x6 uint8 벡터 (S1, S2, S3, S4, S5, S6)
    %       tables:    gf_arithmetic.m에서 생성된 GF 테이블 구조체
    %
    %   출력:
    %       lambda:    오류 위치 다항식 σ(x) (계수 벡터)
    %       omega:     오류 평가 다항식 ω(x) (계수 벡터)

    % --- 1. 초기화 ---
    
    TWO_T = 6; % RS(127, 121)의 2t 값
    
    % 다항식 연산 유틸리티 로드
    poly_utils = gf_poly_utils(tables);
    
    % 알고리즘 상태 변수 초기화 (논문 5.2.3절 초기값 참조 [cite: 1311-1312])
    % 다항식은 [p0, p1, p2, ...] 형태의 계수 벡터로 표현합니다.
    C_poly = uint8(1);      % σ^k(x), 오류 위치 다항식. 초기값 σ^0(x) = 1
    W_poly = uint8(1);      % ω^k(x), 오류 평가 다항식. 초기값 ω^0(x) = 1
    B_poly = uint8(1);      % λ^k(x), 보조 다항식 (σ 백업용). 초기값 λ^0(x) = 1
    Beta_poly = uint8(1);   % β^k(x), 보조 다항식 (ω 백업용). 초기값 β^0(x) = 1
    L = 0;                  % l_k, 유효 길이. 초기값 l_0 = 0
    Gamma = uint8(1);       % γ_k, 스케일링 인자. 초기값 γ_0 = 1
    
    % 입력 신드롬 벡터 (S1 ~ S6)
    S = syndromes;

    % --- 2. BM 알고리즘 2t번 반복 (k=0 부터 5까지) ---
    
    % MATLAB 루프는 1-based 이므로, k_iter = 1...6 (k = 0...5)
    for k_iter = 1:TWO_T
        
        k = k_iter - 1; % 0-based 알고리즘 인덱스 k
        
        % (a) 불일치 값(Discrepancy) δ_{k+1} 계산 [cite: 1312]
        % δ_{k+1} = S_{k+1} + σ_1*S_k + σ_2*S_{k-1} + ... + σ_l*S_{k+1-l}
        Delta = uint8(0);
        
        % L_current는 l_k를 의미 (현재 단계 k에서의 l 값)
        L_current = L; 
        
        for j = 0:L_current
            syndrome_index = (k + 1) - j; % S_{k+1-j}
            sigma_index = j + 1;          % σ_j (벡터는 1-based index)
            
            % 유효한 신드롬 인덱스(S1~S6)와 시그마 계수 범위 내에서만 계산
            if syndrome_index >= 1 && syndrome_index <= length(S) && sigma_index <= length(C_poly)
                S_val = S(syndrome_index);
                C_val = C_poly(sigma_index);
                
                term = poly_utils.gf_multiply(S_val, C_val);
                Delta = poly_utils.gf_add(Delta, term);
            end
        end
        
        % (b) 주요 다항식 σ^{k+1}(x), ω^{k+1}(x) 업데이트 [cite: 1313-1314]
        % σ^{k+1}(x) = γ_k * σ^k(x) - δ_{k+1} * λ^k(x) * x
        % ω^{k+1}(x) = γ_k * ω^k(x) - δ_{k+1} * β^k(x) * x
        
        term1_C = poly_utils.scale(C_poly, Gamma);
        term2_C = poly_utils.scale(poly_utils.shift(B_poly, 1), Delta);
        C_next = poly_utils.add(term1_C, term2_C);
        
        term1_W = poly_utils.scale(W_poly, Gamma);
        term2_W = poly_utils.scale(poly_utils.shift(Beta_poly, 1), Delta);
        W_next = poly_utils.add(term1_W, term2_W);

        % (c) 보조 변수 업데이트 (조건부) [cite: 1315-1327]
        
        if Delta == 0 || (2 * L_current > k)
            % 조건 1: δ_{k+1} = 0  또는 2*l_k > k
            % σ, ω의 구조적인 변경이 필요 없음
            
            B_next = poly_utils.shift(B_poly, 1);       % λ^{k+1}(x) = x * λ^k(x)
            Beta_next = poly_utils.shift(Beta_poly, 1); % β^{k+1}(x) = x * β^k(x)
            L_next = L_current;                         % l_{k+1} = l_k
            Gamma_next = Gamma;                         % γ_{k+1} = γ_k
            
        else
            % 조건 2: δ_{k+1} ≠ 0  그리고 2*l_k ≤ k
            % σ, ω의 구조적인 변경(차수 증가)이 필요함
            
            B_next = C_poly;                            % λ^{k+1}(x) = σ^k(x)
            Beta_next = W_poly;                         % β^{k+1}(x) = ω^k(x)
            L_next = (k + 1) - L_current;               % l_{k+1} = k + 1 - l_k
            Gamma_next = Delta;                         % γ_{k+1} = δ_{k+1}
        end
        
        % (d) 다음 반복을 위해 상태 저장
        C_poly = C_next;
        W_poly = W_next;
        B_poly = B_next;
        Beta_poly = Beta_next;
        L = L_next;
        Gamma = Gamma_next;
        
    end % end of k_iter loop

    % --- 3. 최종 결과 반환 ---
    lambda = C_poly;
    omega = W_poly;

end
