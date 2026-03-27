clear
rng(0)

% ---------------- Fixed parameters ----------------
EbNo_dB_fixed = 10.5;     % choose 10 or 10.5
T_list = 15:2:30;           % sweep T from 15 to 30

eta = 1;
w   = 25;
flip_num = 3;

p = 4;                    % bits per symbol
q = 2^p;

load('arith_16.mat');
load('204.102.3.6.16.mat');
load('G_204_102.mat');

h = full(h);
N = size(h,2);
M = size(h,1);
K = N - M;
R = K/N;

CN_lst = cell(M,1);
for i = 1:M
    CN_lst{i,1} = find(h(i,:));
end

info_seq = randi([0 q-1], 1, K);
code_seq = gf_mat_mul(info_seq, G, add_mat, mul_mat);

% ---------------- QAM mapping ----------------

qam16 = [
    -3 + 3i;  % 0011
    -3 + 1i;  % 0010
    -3 - 1i;  % 0001
    -3 - 3i;  % 0000
    -1 + 3i;  % 0111
    -1 + 1i;  % 0110
    -1 - 1i;  % 0101
    -1 - 3i;  % 0100
    1 + 3i;  % 1011
    1 + 1i;  % 1010
    1 - 1i;  % 1001
    1 - 3i;  % 1000
    3 + 3i;  % 1111
    3 + 1i;  % 1110
    3 - 1i;  % 1101
    3 - 3i   % 1100
    ];

qam_binary_map = [ 0  0  0  0;   % -3-3j
                   0  0  0  1;   % -3-1j
                   0  0  1  1;   % -3+1j
                   0  0  1  0;   % -3+3j
                   0  1  1  0;   % -1+3j
                   0  1  1  1;   % -1+1j
                   0  1  0  1;   % -1-1j
                   0  1  0  0;   % -1-3j
                   1  1  0  0;   % +1-3j
                   1  1  0  1;   % +1-1j
                   1  1  1  1;   % +1+1j
                   1  1  1  0;   % +1+3j
                   1  0  1  0;   % +3+3j
                   1  0  1  1;   % +3+1j
                   1  0  0  1;   % +3-1j
                   1  0  0  0];  % +3-3j


avg_pow = qam16' * qam16 / q;
nrm_fct = sqrt(avg_pow);
gf16 = 0:q-1;
%alph_bin =  logical(fliplr(dec2bin(gf16, p) - 48)); % symbols in binary

% ---------------- Storage for T sweep ----------------
numT = length(T_list);

FE         = zeros(numT,1);
genFrame   = zeros(numT,1);
iters_cnr  = zeros(numT,1);
BE_Coded   = zeros(numT,1);
BE_unCoded = zeros(numT,1);

FER_vs_T   = zeros(numT,1);
BER_vs_T   = zeros(numT,1);
AvgIt_vs_T = zeros(numT,1);

targetFE = 500;
max_gen   = 1e5;

% ---------------- Fixed SNR quantities ----------------
EbNo_linear = 10^(EbNo_dB_fixed/10);
avg_symbol_energy = 1;
No = avg_symbol_energy / (p * EbNo_linear * R);
sigma0 = sqrt(No/2) * nrm_fct;
nse_std = eta * sigma0;

% ---------------- Sweep over T ----------------
for t_idx = 1:numT

    T = T_list(t_idx);

    while (FE(t_idx) < targetFE && genFrame(t_idx) < max_gen)

        genFrame(t_idx) = genFrame(t_idx) + 1;

        c = qam16(code_seq' + 1).';
        

        n = sigma0*randn(1,N) + 1i*sigma0*randn(1,N);
        y = c + n;

        hard_d_cmplx = zeros(1,N);
        hard_d_gf16  = zeros(1,N);

        for j = 1:N
            distance = abs(qam16 - y(j));
            [~, min_idx] = min(distance);
            hard_d_cmplx(j) = qam16(min_idx);
            hard_d_gf16(j)  = gf16(min_idx);
        end

        % -------- uncoded BER --------

         errors_uncoded_bit = zeros(1, K);

        for e = 1 : K          
                if hard_d_gf16(e)~=code_seq(e) 
                    s1 = dec2bin(code_seq(e),p);
                    s2 = dec2bin(hard_d_gf16(e),p);
                    code_seq_ = double(s1);
                    dec_seq_ = double(s2);
                    num_diff_bit = sum(code_seq_ ~= dec_seq_);
                    errors_uncoded_bit(e) = errors_uncoded_bit(e) + num_diff_bit;                
                end
        end
        
        un_bit_error = sum(errors_uncoded_bit); % no of bit error in each frame
        BE_unCoded(t_idx) = BE_unCoded(t_idx) + un_bit_error;

        % -------- decoder --------
        [seqgf, failed, l] = decodeMultivote(code_seq, hard_d_cmplx, hard_d_gf16, ...
            qam16, gf16, y, h, N, M, T, w, add_mat, mul_mat, div_mat, ...
            CN_lst, nse_std, qam_binary_map, flip_num, No);

        % -------- coded BER --------
         errors_coded_bit = zeros(1, K);
            for g = 1 : K            
                if seqgf(g)~=code_seq(g)
                    s1 = dec2bin(code_seq(g),p);
                    s2 = dec2bin(seqgf(g),p);
                    code_seq_ = double(s1);
                    dec_seq_ = double(s2);
                    num_diff_bit = sum(code_seq_ ~= dec_seq_);
                    errors_coded_bit(g) = errors_coded_bit(g) + num_diff_bit;                
                end     
            end

        bit_error = sum(errors_coded_bit);
        BE_Coded(t_idx) = BE_Coded(t_idx) + bit_error;

        iters_cnr(t_idx) = iters_cnr(t_idx) + l;

        if bit_error > 0
            FE(t_idx) = FE(t_idx) + 1;
        end

    end

    FER_vs_T(t_idx)   = FE(t_idx) / genFrame(t_idx);
    BER_vs_T(t_idx)   = BE_Coded(t_idx) / (genFrame(t_idx) * K * p);
    AvgIt_vs_T(t_idx) = iters_cnr(t_idx) / genFrame(t_idx);

    fprintf('T=%d, Eb/No=%.1f dB: FER=%.3e, BER=%.3e, AvgIt=%.2f\n', ...
        T, EbNo_dB_fixed, FER_vs_T(t_idx), BER_vs_T(t_idx), AvgIt_vs_T(t_idx));
end

figure;
semilogy(T_list, BER_vs_T, 'o-', 'LineWidth', 1.2);
grid on;
xlabel('T');
ylabel('BER');
title(sprintf('BER vs T at E_b/N_0 = %.1f dB', EbNo_dB_fixed));

figure;
semilogy(T_list, FER_vs_T, 's-', 'LineWidth', 1.2);
grid on;
xlabel('T');
ylabel('FER');
title(sprintf('FER vs T at E_b/N_0 = %.1f dB', EbNo_dB_fixed));


figure;
plot(T_list, AvgIt_vs_T, 'd-', 'LineWidth', 1.2);
grid on;
xlabel('T');
ylabel('Average Iterations');
title(sprintf('Average Iterations vs T at E_b/N_0 = %.1f dB', EbNo_dB_fixed));

