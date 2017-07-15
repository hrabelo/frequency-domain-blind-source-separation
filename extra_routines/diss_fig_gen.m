function diss_fig_gen(fig_num, varargin)
% Função para gerar figuras para a dissertação
%
% diss_fig_gen(1) - SUBGAUSSIAN VS SUPERGAUSSIAN
% diss_fig_gen(2) - BSS AMBIGUITIES
% diss_fig_gen(3) - SCATTER ICA
% diss_fig_gen(4) - SOURCE DIST
% diss_fig_gen(5) - SPEECH SPECTROGRAM
% diss_fig_gen(6) - SPEECH FREQUENCY ENVELOPE
% diss_fig_gen(7) - CORRELATION SPEAKER
% diss_fig_gen(8) - WINDOW SPECTRUM
% diss_fig_gen(9, W, d, fs) - DIRECTIVITY PATTERN
% diss_fig_gen(10, W, d, fs) - DIRECTIVITY PATTERN 3 SOURCES
% diss_fig_gen(12, W, d, fs) - DIRECTIVITY PATTERN MEAN
% diss_fig_gen(13, W, d, fs) - DOA
% diss_fig_gen(14, W, Y, J, fs) - POW SPECTROGRAM
% diss_fig_gen(15, W, Y, fs) - POW CORRELATION
% diss_fig_gen(16) - GEN GAUSSIAN
% diss_fig_gen(17) - GEN GAUSSIAN R
% diss_fig_gen(18) - COMPARISON PLOTS
% diss_fig_gen(19) - REVERB DECAY
% diss_fig_gen(20) - SCATTER ICA 2
%

switch fig_num
    
    case 1
        %% SUBGAUSSIAN VS SUPERGAUSSIAN
        x = -5:0.01:5;
        % Laplace
        f = 0.5 * exp(-abs(x));
        plot(x, f, '--k', 'LineWidth', 1);
        xlabel('x'); ylabel('q(x)'); ylim([0 0.55]);
        hold;
        
        % Gaussiana
        f = 1/sqrt(2*pi) * exp(-0.5*x.^2);
        plot(x, f, '-k', 'LineWidth', 1);
        
        % Wigner Semicircle
        f = 1/(2*pi) * sqrt(4 - x.^2);
        plot(x, f, ':k', 'LineWidth', 1.5);

        legend('Supergaussiana', 'Gaussiana', 'Subgaussiana');
        

    case 2
        %% BSS AMBIGUITIES
        secs = 10;
        ts = 0.01;
        T = 1.5;
        factor = 1.5;
        scaling = [0.8 0.2 -1.2];
        perm = [1 3 2];
        x_range = [0 5.5*T/ts];
        
        s(1,:) = cos((2*pi/T*(0:ts:secs)));
        s(2,:) = square((2*pi/T*(0:ts:secs)));
        s(3,:) = sawtooth((2*pi/T*(0:ts:secs)));
        
        h = [0.5 1.2 1.1; 1.2 0.3 0.6; 1.5 0.1 0.9];
        x = h * s;
        w = diag(scaling)*inv(h);
        w = w(perm, :);
        y = w*x;
        
        subplot(3,3,1); plot(s(1,:), 'k', 'LineWidth', 1);
        y_range = factor*[min(s(1,:)) max(s(1,:))];
        set(gca, 'XTick', [], 'YTick', []); xlim(x_range); ylim(y_range);
        xlabel('Fonte s_1');

        subplot(3,3,4); plot(s(2,:), 'k', 'LineWidth', 1); 
        set(gca, 'XTick', [], 'YTick', []); xlim(x_range); ylim(y_range);
        xlabel('Fonte s_2');

        subplot(3,3,7); plot(s(3,:), 'k', 'LineWidth', 1); 
        set(gca, 'XTick', [], 'YTick', []); xlim(x_range); ylim(y_range);
        xlabel('Fonte s_3');

        subplot(3,3,3); plot(y(1,:), 'k', 'LineWidth', 1);
        set(gca, 'XTick', [], 'YTick', []); xlim(x_range); ylim(y_range);
        xlabel('Saída y_1');

        subplot(3,3,6); plot(y(2,:), 'k', 'LineWidth', 1);
        set(gca, 'XTick', [], 'YTick', []); xlim(x_range); ylim(y_range);
        xlabel('Saída y_2');

        subplot(3,3,9); plot(y(3,:), 'k', 'LineWidth', 1);
        set(gca, 'XTick', [], 'YTick', []); xlim(x_range); ylim(y_range);
        xlabel('Saída y_3');

        subplot(3,3,2); plot(x(1,:), 'k', 'LineWidth', 1); 
        y_range = factor*[min(x(2,:)) max(x(2,:))];
        set(gca, 'XTick', [], 'YTick', []); xlim(x_range); ylim(y_range);
        xlabel('Mistura x_1');

        subplot(3,3,5); plot(x(2,:), 'k', 'LineWidth', 1);
        set(gca, 'XTick', [], 'YTick', []); xlim(x_range); ylim(y_range);
        xlabel('Mistura x_2');

        subplot(3,3,8); plot(x(3,:), 'k', 'LineWidth', 1);
        set(gca, 'XTick', [], 'YTick', []); xlim(x_range); ylim(y_range);
        xlabel('Mistura x_3');
        
    case 3
        %% SCATTER ICA 
        
        Path = 'C:\Users\Victorio\Documents\MATLAB\dev\';
        num_samp = 32000;
        
        s(:,1) = wavread([Path 'female_src_2.wav'], num_samp); 
        s(:,2) = wavread([Path 'male_src_2.wav'], num_samp);
        s = s.';
        
        h = rand(2, 2);
        x = h * s;
        
        [z v] = pre_whitening(x);
        y = natICA(x, 'InitSepMat', v, 'InitSourceSig', z, 'ScoreFunction', 'sign');
        
        figure; scatter(s(1,:), s(2,:), 2, 'k');
        set(gca, 'XTick', [], 'YTick', [], 'color', 'white'); xlabel('s_1'); ylabel('s_2');
        
        figure; scatter(x(1,:), x(2,:), 2, 'k');
        set(gca, 'XTick', [], 'YTick', [], 'color', 'white'); xlabel('x_1'); ylabel('x_2');
        
        figure; scatter(z(1,:), z(2,:), 2, 'k');
        set(gca, 'XTick', [], 'YTick', [], 'color', 'white'); xlabel('z_1'); ylabel('z_2');
        
        figure; scatter(y(1,:), y(2,:), 2, 'k');
        set(gca, 'XTick', [], 'YTick', [], 'color', 'white'); xlabel('y_1'); ylabel('y_2');
        
    case 4
        %% SOURCE DISTS
        
        Path = 'C:\Users\Victorio\Documents\MATLAB\dev\';
        num_samp = 96000;
       
        s(:,1) = wavread([Path 'female_src_1.wav'], num_samp);  s(:,1) = (s(:,1) - mean(s(:,1))) ./ std(s(:,1));
        s(:,2) = wavread([Path 'female_src_3.wav'], num_samp);  s(:,2) = (s(:,2) - mean(s(:,2))) ./ std(s(:,2));
        s(:,3) = wavread([Path 'male_src_1.wav'], num_samp);    s(:,3) = (s(:,3) - mean(s(:,3))) ./ std(s(:,3));
        s(:,4) = wavread([Path 'male_src_2.wav'], num_samp);    s(:,4) = (s(:,4) - mean(s(:,4))) ./ std(s(:,4));
        s = resample(s,1,2);
        s = s.';
        
        x(1,:) = s(1,:) + s(3,:);                               x(1,:) = (x(1,:) - mean(x(1,:))) ./ std(x(1,:));
        x(2,:) = s(1,:) + s(2,:) + s(3,:) + s(4,:);             x(2,:) = (x(2,:) - mean(x(2,:))) ./ std(x(2,:));
         
        [b,a] = hist(s(1,:), 250); b = b/trapz(a,b);
        g = 1/sqrt(2*pi) * exp(-0.5*a.^2);
        figure;
        plot(a, b, 'k', 'LineWidth', 1); hold; plot(a, g, '--k', 'LineWidth', 1);
        xlabel('Fonte 1', 'FontSize', 14); xlim([-4 4]);
    
        [b,a] = hist(x(1,:), 250); b = b/trapz(a,b);
        g = 1/sqrt(2*pi) * exp(-0.5*a.^2);
        figure;
        plot(a, b, 'k', 'LineWidth', 1); hold; plot(a, g, '--k', 'LineWidth', 1);
        xlabel('Fonte 1 + Fonte 3', 'FontSize', 14);  xlim([-5 5]);

        [b,a] = hist(x(2,:), 250); b = b/trapz(a,b);
        g = 1/sqrt(2*pi) * exp(-0.5*a.^2);
        figure;
        plot(a, b, 'k', 'LineWidth', 1); hold; plot(a, g, '--k', 'LineWidth', 1);
        xlabel('Fonte 1 + Fonte 2 + Fonte 3 + Fonte 4', 'FontSize', 14); xlim([-7 7]);

    case 5
        %% SPEECH SPECTROGRAM
        Path = 'C:\Documents and Settings\Victorio\Meus documentos\MATLAB\dev\';
        num_samp = 96000;
       
        [s fs] = wavread([Path 'female_src_1.wav'], num_samp);
        s = resample(s,1,2);    fs = fs/2;
        s = s.';
        
        K = 1024;
        windlen = 512;
        J = windlen/4;
        wind = 0.5*hann(windlen, 'periodic').';
        S = squeeze(stft(s, K, J, wind, 'zeropad'));
        
        freq = [0 (fs/2 + fs/K)];
        time = [0 (size(S, 2) - 1)*(J/fs)];
        time2 = (0:size(s, 2) - 1)*(1/fs);
        
        envS = envelope(S(1:K/2 + 1, :), 4);

        figure; plot(time2, s, 'k');
        set(gca, 'YTick', []); xlabel('Tempo (segundos)');
        
        figure;        colormap('gray');
        imagesc(time, freq, envS);  set(gca, 'Ydir', 'normal');
        xlabel('Tempo (segundos)'); ylabel('Frequência (Hertz)');
        
    case 6
        %% SPEECH FREQUENCY ENVELOPE
        Path = 'C:\Documents and Settings\Victorio\Meus documentos\MATLAB\dev\';
        num_samp = 96000;
       
        [s fs] = wavread([Path 'female_src_1.wav'], num_samp);
        s = resample(s,1,2);    fs = fs/2;
        s = s.';
        
        K = 1024;
        windlen = 512;
        J = windlen/4;
        wind = 0.5*hann(windlen, 'periodic').';
        S = squeeze(stft(s, K, J, wind, 'zeropad'));
        
        f1 = 56;
        f2 = 57;
        f3 = 82;
        f4 = 163;
        
        freq = 0:fs/K:fs/2 + fs/K;
        time = (0:size(S, 2) - 1)*(J/fs);
        
        envS = envelope(S(1:K/2 + 1, :), 1);

        subplot(2,2,1); plot(time, envS(f1, :), 'k', 'LineWidth', 1);
        set(gca, 'YTick', [], 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        xlabel('Tempo (segundos)'); title([num2str(freq(f1), 4) ' Hz'], 'FontSize', 14);
        
        subplot(2,2,2); plot(time, envS(f2, :), 'k', 'LineWidth', 1);
        set(gca, 'YTick', [], 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        xlabel('Tempo (segundos)'); title([num2str(freq(f2), 4) ' Hz'], 'FontSize', 14);
        
        subplot(2,2,3); plot(time, envS(f3, :), 'k', 'LineWidth', 1);
        set(gca, 'YTick', [], 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        xlabel('Tempo (segundos)'); title([num2str(freq(f3), 4) ' Hz'], 'FontSize', 14);
        
        subplot(2,2,4); plot(time, envS(f4, :), 'k', 'LineWidth', 1);
        set(gca, 'YTick', [], 'Position', get(gca, 'OuterPosition') - ...
            get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
        xlabel('Tempo (segundos)'); title([num2str(freq(f4), 4) ' Hz'], 'FontSize', 14);
        
    case 7
        %% CORRELATION SPEAKER
        Path = 'C:\Documents and Settings\Victorio\Meus documentos\MATLAB\dev\';
        num_samp = 96000;
       
        [s(:, 1) fs] = wavread([Path 'female_src_1.wav'], num_samp);
        s(:, 2) = wavread([Path 'male_src_3.wav'], num_samp);
        s = resample(s,1,2);    fs = fs/2;
        s = s.';

        K = 4096;
        windlen = 2048;
        J = windlen/4;
        wind = 0.5*hann(windlen, 'periodic').';
        S = stft(s, K, J, wind, 'zeropad');

        freq = [0 (fs/2 + fs/K)];
        
        tmp1 = envelope(squeeze(S(1, 1:K/2+1, :)), 1);
        mR1 = mean(tmp1, 2);
        vR1 = var(tmp1, 0, 2);

        tmp2 = envelope(squeeze(S(2, 1:K/2+1, :)), 1);
        mR2 = mean(tmp2, 2);
        vR2 = var(tmp2, 0, 2);

        RR = fast_corr2by2(tmp1, tmp1, mR1, mR1, vR1, vR1);
        figure;
        colormap('gray');
        imagesc(freq, freq, RR, [-0.4 1]);  set(gca, 'Ydir', 'normal');
        xlabel('Frequência (Hertz)'); ylabel('Frequência (Hertz)');
        tmptmp = RR(100:112, 100:112);
        disp(tmptmp)
        
        RR = fast_corr2by2(tmp1, tmp2, mR1, mR2, vR1, vR2);
        figure;
        colormap('gray');
        imagesc(freq, freq, RR, [-0.4 1]);  set(gca, 'Ydir', 'normal');
        xlabel('Frequência (Hertz)'); ylabel('Frequência (Hertz)');
        tmptmp = RR(100:112, 100:112);
        disp(tmptmp)

    case 8
        %% WINDOW SPECTRUM
        windlen = 1024;
        totallen = 2^21;
        lim_Lf = 10;
        lim_db = -150;
        
        
        samp_per_bin = totallen/windlen;
        Lf = (0:1/samp_per_bin:windlen - 1/samp_per_bin);
        maxLf_ind = find(Lf > lim_Lf, 1);
        Lf = Lf(1:maxLf_ind);

        wind = ones(1,windlen);
        windf = (1/length(wind))*fft(wind, totallen);
        windf = 10*log10(windf.*conj(windf));
        windf = windf - max(windf);
        windf = windf(1:maxLf_ind);
        
        figure;
        plot(Lf, windf, 'k', 'LineWidth', 1); ylim([lim_db 2]); xlim([0 lim_Lf]);
        xlabel('$$\frac{f}{f_s/L}$$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$$10 \times \log(|win|^{2})$$', 'Interpreter', 'latex', 'FontSize', 12);
        
        wind = 0.5*hann(windlen, 'periodic').';
        windf = (1/length(wind))*fft(wind, totallen);
        windf = 10*log10(windf.*conj(windf));
        windf = windf - max(windf);
        windf = windf(1:maxLf_ind);
        
        figure;
        plot(Lf, windf, 'k', 'LineWidth', 1); ylim([lim_db 2]); xlim([0 lim_Lf]);
        xlabel('$$\frac{f}{f_s/L}$$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$$10 \times \log(|win|^{2})$$', 'Interpreter', 'latex', 'FontSize', 12);
        
        wind = 0.72*window(@blackmanharris, windlen).';
        windf = (1/length(wind))*fft(wind, totallen);
        windf = 10*log10(windf.*conj(windf));
        windf = windf - max(windf);
        windf = windf(1:maxLf_ind);
        
        figure;
        plot(Lf, windf, 'k', 'LineWidth', 1); ylim([lim_db 2]); xlim([0 lim_Lf]);
        xlabel('$$\frac{f}{f_s/L}$$', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$$10 \times \log(|win|^{2})$$', 'Interpreter', 'latex', 'FontSize', 12);
        
        figure;     sidelen = floor(0.1*windlen);
        wind = [zeros(1, sidelen)   hann(windlen, 'periodic').'         zeros(1, sidelen)]; plot(wind, 'k', 'LineWidth', 1); hold;
        wind = [zeros(1, sidelen)   window(@blackmanharris, windlen).'  zeros(1, sidelen)]; plot(wind, '--k', 'LineWidth', 1);
        wind = [zeros(1, sidelen)   kaiser(windlen, 18).'               zeros(1, sidelen)]; plot(wind, ':k', 'LineWidth', 1);
        wind = [zeros(1, sidelen)   ones(1, windlen)                    zeros(1, sidelen)]; plot(wind, 'k', 'LineWidth', 1);
        ylim([0 1.3]); xlim([0 length(wind)]);
        set(gca, 'XTick', []);
        legend('Hanning', 'Blackman-Harris', 'Kaiser-Bessel (\beta = 18)');


        
    case 9
        %% DIRECTIVITY PATTERN
        W = varargin{1};
        d = varargin{2};
        fs = varargin{3};
        K = size(W, 1);
        M = size(W, 3);
        num_teta = 360;
        
        freq = [0 : fs/K : (fs/K)*(K/2-1), -fs/2 : fs/K : -fs/K];
        pos = -(M-1)*d/2 : d : (M-1)*d/2;
        %pos = 0: d : (M-1)*d;
        teta = 0:pi/num_teta:pi;
        tetad = teta*180/pi;

        for i = [round(0.08*K) round(0.18*K) round(0.4*K)]
            U = dirpattern(squeeze(W(i, :, :)), pos, freq(i), teta);
            
            figure; plot(tetad, U(1,:), 'k', 'LineWidth', 1); hold;
            plot(tetad, U(2,:), '--k', 'LineWidth', 1);
            xlabel('$$\theta$$ (em graus)', 'Interpreter', 'latex', 'FontSize', 14);
            title([num2str(freq(i), 4) ' Hz'], 'FontSize', 16);
            set(gca, 'YTick', []);         legend('Fonte 1', 'Fonte 2');
        end
        
    case 10
        %% DIRECTIVITY PATTERN 3 SOURCES
        W = varargin{1};
        d = varargin{2};
        fs = varargin{3};
        K = size(W, 1);
        M = size(W, 3);
        num_teta = 360;
        
        freq = [0 : fs/K : (fs/K)*(K/2-1), -fs/2 : fs/K : -fs/K];
        pos = -(M-1)*d/2 : d : (M-1)*d/2;
        teta = 0:pi/num_teta:pi;
        tetad = teta*180/pi;

        i = round(0.21*K);
        U = dirpattern(squeeze(W(i, :, :)), pos, freq(i), teta);

        figure; plot(tetad, U(1,:), 'k', 'LineWidth', 1); hold;
        plot(tetad, U(2,:), '--k', 'LineWidth', 1);
        plot(tetad, U(3,:), ':k', 'LineWidth', 1);
        xlabel('$$\theta$$ (em graus)', 'Interpreter', 'latex', 'FontSize', 14);
        title([num2str(freq(i), 4) ' Hz'], 'FontSize', 16);
        set(gca, 'YTick', []);         legend('Fonte 1', 'Fonte 2', 'Fonte 3');
        
    case 11
        %% DIRECTIVITY PATTERN MEAN
        W = varargin{1};
        d = varargin{2};
        fs = varargin{3};
        K = size(W, 1);
        N = size(W, 2);
        M = size(W, 3);
        num_teta = 360;
        
        freq = [0 : fs/K : (fs/K)*(K/2-1), -fs/2 : fs/K : -fs/K];
        pos = -(M-1)*d/2 : d : (M-1)*d/2;
        teta = 0:pi/num_teta:pi;
        tetad = teta*180/pi;

        U = zeros(N, length(teta));
        for i = 1:K
            Utmp = dirpattern(squeeze(W(i, :, :)), pos, freq(i), teta);
            U = U + Utmp;
        end
        U = U / (K/2 + 1);
        
        figure; plot(tetad, U(1,:), 'k', 'LineWidth', 1); hold;
        plot(tetad, U(2,:), '--k', 'LineWidth', 1);
        plot(tetad, U(3,:), ':k', 'LineWidth', 1);
        xlabel('$$\theta$$ (em graus)', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca, 'YTick', []);         legend('Fonte 1', 'Fonte 2', 'Fonte 3');
        
    case 12
        %% DOA
        W = varargin{1};
        d = varargin{2};
        fs = varargin{3};
        K = size(W, 1);
        N = size(W, 2);
         
        freq = 0 : fs/K : (fs/K)*(K/2);

        teta = zeros(N, K/2+1);
        for i = 1:K/2+1
            teta(:, i) = doa(squeeze(W(i, :, :)), d, freq(i));
        end
        
        teta = teta*180/pi;
        nanteta = (sum(isnan(teta), 1) > 0);
        disp(sum(nanteta));
        mteta = median(teta(:, ~nanteta), 2);
        
        figure; plot(freq, teta(1, :), 'xk'); hold;
        plot(freq, teta(2, :), '*b');
        if N > 2
            plot(freq, teta(3, :), '+r');
            if N > 3
                plot(freq, teta(4, :), 'og');
            end
        end
        plot(freq, mteta(1)*ones(1,K/2+1), 'k', 'LineWidth', 1);
        plot(freq, mteta(2)*ones(1,K/2+1), 'b', 'LineWidth', 1);
        if N > 2
            plot(freq, mteta(3)*ones(1,K/2+1), 'r', 'LineWidth', 1);
            if N > 3
                plot(freq, mteta(4)*ones(1,K/2+1), 'g', 'LineWidth', 1);
            end
        end
        xlabel('Frequência (em Hertz)', 'FontSize', 12);    ylim([0 180]);
        ylabel('$$\theta$$ (em graus)', 'Interpreter', 'latex', 'FontSize', 14);
        if N == 2,      legend('Fonte 1', 'Fonte 2');
        elseif N == 3,  legend('Fonte 1', 'Fonte 2', 'Fonte 3');
        elseif N == 4,  legend('Fonte 1', 'Fonte 2', 'Fonte 3', 'Fonte 4');
        end
        disp(mteta);

    case 14
        %% POW SPECTROGRAM
        W = varargin{1};
        Y = varargin{2};
        J = varargin{3};
        fs = varargin{4};
        
        K = size(W,1);
        
        freq = [0 (fs/2 + fs/K)];
        time = [0 (size(Y, 2) - 1)*(J/fs)];
        
        envY = zeros(size(Y, 1), K/2 + 1, size(Y, 3));
        for ck = 1:K/2 + 1
            envY(:, ck, :) = envelope(squeeze(Y(:, ck, :)), 3, squeeze(W(ck, :, :)));
        end
        
        figure;        colormap('gray');
        imagesc(time, freq, squeeze(envY(1, :, :)), [0 1]);  set(gca, 'Ydir', 'normal');
        xlabel('Tempo (segundos)'); ylabel('Frequência (Hertz)'); title('powRatio_1');

        figure;        colormap('gray');
        imagesc(time, freq, squeeze(envY(2, :, :)), [0 1]);  set(gca, 'Ydir', 'normal');
        xlabel('Tempo (segundos)'); ylabel('Frequência (Hertz)'); title('powRatio_2');

    case 15
        %% POW CORRELATION
        W = varargin{1};
        Y = varargin{2};
        fs = varargin{3};
        
        K = size(W,1);
        
        freq = [0 (fs/2 + fs/K)];

        envY = zeros(size(Y, 1), K/2 + 1, size(Y, 3));
        for ck = 1:K/2 + 1
            envY(:, ck, :) = envelope(squeeze(Y(:, ck, :)), 3, squeeze(W(ck, :, :)));
        end
        
        tmp1 = squeeze(envY(1, :, :));
        mR1 = mean(tmp1, 2);
        vR1 = var(tmp1, 0, 2);

        tmp2 = squeeze(envY(2, :, :));
        mR2 = mean(tmp2, 2);
        vR2 = var(tmp2, 0, 2);

        RR = fast_corr2by2(tmp1, tmp1, mR1, mR1, vR1, vR1);
        figure;
        colormap('gray');
        imagesc(freq, freq, RR, [-0.4 1]);  set(gca, 'Ydir', 'normal');
        xlabel('Frequência (Hertz)'); ylabel('Frequência (Hertz)');
        
        RR = fast_corr2by2(tmp1, tmp2, mR1, mR2, vR1, vR2);
        figure;
        colormap('gray');
        imagesc(freq, freq, RR, [-0.4 1]);  set(gca, 'Ydir', 'normal');
        xlabel('Frequência (Hertz)'); ylabel('Frequência (Hertz)');
        
    case 16
        %% GEN GAUSSIAN
        range = -4:0.01:4;
        lr = length(range);
        y = zeros(lr, lr);
        
        for cl = 1:lr
            y(cl, :) = complex(range, ones(1, lr)*range(cl));
        end


        sigma = 1;
        r = 0.5;
        genG = (r / (2*sigma*gamma(1/r)))*exp(-(1/r)*abs(y / sigma).^r);
        figure;   mesh(range, range, genG); colormap jet
        xlabel('Parte real'); ylabel('Parte imaginária');

        r = 1;
        genG = (r / (2*sigma*gamma(1/r)))*exp(-(1/r)*abs(y / sigma).^r);
        figure;   mesh(range, range, genG); colormap jet
        xlabel('Parte real'); ylabel('Parte imaginária');

        r = 4;
        genG = (r / (2*sigma*gamma(1/r)))*exp(-(1/r)*abs(y / sigma).^r);
        figure;   mesh(range, range, genG); colormap jet
        xlabel('Parte real'); ylabel('Parte imaginária');

    case 17
        %% GEN GAUSSIAN R
        r = 0.2:0.001:2;
        kurt = gamma(5./r) .* gamma(1./r) ./ (gamma(3./r)).^2;
        figure;
        semilogy(r, kurt, 'k', 'LineWidth', 1);
        xlabel('r', 'FontSize', 12); ylabel('Curtose', 'FontSize', 12); 
        
        r = 2:0.001:14;
        kurt = gamma(5./r) .* gamma(1./r) ./ (gamma(3./r)).^2;
        figure;
        plot(r, kurt, 'k', 'LineWidth', 1);
        xlabel('r', 'FontSize', 12); ylabel('Curtose', 'FontSize', 12); 
  
    case 18
        %% COMPARISON PLOTS

        % CLUSTER2D
        sup = [32.3, 32.5, 30.1, 28.2, 24.1, 20.5, 16];
        tdoa = [32.3, 32.5, 27.2, 19.9, 15.7, 10.8, 9.1];
        doaconj = [32.3, 32.3, 22, 16.8, 13.1, 10.6, 6.6];
        doaharm = [32.3, 32.3, 22, 16.6, 13, 8.6, 6.4];
        doaall = [31.9, 30, 18.8, 14.9, 12.1, 10.6, 7];
        conjcorr = [2 2.4, 2.4, 3.7, 3.2, 3, 2.5];
        all = [30.5, 29.2, 26.4, 20.8, 20.6, 19.6, 15.2];
        t60 = [0 50 100 150 200 250 500];

        % LINE2D ERRADO (ângulos iguais ao cluster)
        sup_e = [31.9, 30.6, 26.6, 25.6, 22.9, 21.2, 17.7];
        tdoa_e = [20.7, 25.3, 21.1, 19, 15.8, 13.2, 9.1];
        doaconj_e = [31.9, 29.8, 21, 18.6, 17.4, 14.2, 10];
        doaharm_e = [31.9, 29.7, 20.7, 17.3, 15.8, 14, 9.1];
        doaall_e = [31.9, 30.1, 24.7, 20, 17.9, 15.3, 10.7];
        conjcorr_e = [3.4 3.7, 3.9, 4, 4, 4.4, 2.9];
        all_e = [31.9, 30.3, 25.3, 25.1, 22.6, 20.8, 15.7];
        
        % LINE2D
        sup_l = [32.9, 32.4, 30.4, 26.9, 22, 19.9, 15.7];
        tdoa_l = [29.2, 28.9, 26.2, 16.3, 9.9, 8.8, 5.2];
        doaconj_l = [32.9, 32.4, 30.3, 16.4, 10, 8.1, 7.1];
        doaharm_l = [32.9, 32.4, 24.1, 16.4, 9.4, 7.5, 6.5];
        doaall_l = [32.9, 32.4, 27.4, 16.7, 11, 8.3, 7.5];
        conjcorr_l = [5.1 3.2, 2.8, 3.8, 5.4, 3.5, 2.9];
        all_l = [32.4, 32.4, 27.4, 25.1, 20.6, 19.1, 13.2];
        
        figure;
        plot(t60, sup_l, 'ks-.', 'LineWidth', 1); hold;
        plot(t60, doaconj_l, 'k+:', 'LineWidth', 1);
        plot(t60, doaharm_l, 'kx--', 'LineWidth', 1);
        plot(t60, doaall_l, 'ko-', 'LineWidth', 1);
        ylim([0 34]); legend('Supervisionado', 'DOA + ConjCorr', 'DOA + HarmCorr', 'DOA + GlobalCorr + LocalCorr');
        ylabel('SIR', 'FontSize', 12); xlabel('T_{60} (em milissegundos)', 'FontSize', 12);
        set(gca, 'XTick', [0 50 100 150 200 250 500]);

        figure;
        plot(t60, sup, 'ks-.', 'LineWidth', 1); hold;
        plot(t60, conjcorr, 'kx-', 'LineWidth', 1);
        ylim([0 34]); legend('Supervisionado', 'ConjCorr');
        ylabel('SIR', 'FontSize', 12); xlabel('T_{60} (em milissegundos)', 'FontSize', 12);
        set(gca, 'XTick', [0 50 100 150 200 250 500]);     

        figure;
        plot(t60, doaall, 'ks:', 'LineWidth', 1); hold;
        plot(t60, doaall_e, 'kx-', 'LineWidth', 1);
        ylim([0 34]); legend('Arranjo em cluster', 'Arranjo em linha');
        ylabel('SIR', 'FontSize', 12); xlabel('T_{60} (em milissegundos)', 'FontSize', 12);
        set(gca, 'XTick', [0 50 100 150 200 250 500]);   
        
        figure;
        plot(t60, sup, 'ks-.', 'LineWidth', 1); hold;
        plot(t60, tdoa, 'k+:', 'LineWidth', 1);
        plot(t60, doaall_e, 'kx--', 'LineWidth', 1);
        plot(t60, all, 'ko-', 'LineWidth', 1);
        ylim([0 34]); legend('Supervisionado', 'TDOAclust', 'DOA + GlobalCorr + LocalCorr', 'GlobalCorr + LocalCorr');
        ylabel('SIR', 'FontSize', 12); xlabel('T_{60} (em milissegundos)', 'FontSize', 12);
        set(gca, 'XTick', [0 50 100 150 200 250 500]);

    case 19
        %% REVERB DECAY
        
        mixing_mode = 'ism';
        N = 1;
        M = 1;
        num_samp = 48000;
        reverb_time = 0.8;
        mics_struct = 'cluster2d';
        dist_scrmic = 1;
        ang_src(1) = 90;
        Source{1} = 'C:\Program Files\MATLAB71\work\revtime.wav';
        bss_read;
        t = (0:1/fs:length(s)/fs-1/fs)*1000;
        
        figure; plot(t, x, 'k', 'LineWidth', 1); xlim([0 1000]); ylim([-0.03 0.03]);
        set(gca, 'YTick', []); xlabel('Tempo em milissegundos', 'FontSize', 12);

    case 20
        %% SCATTER ICA 2
        
        Path = 'C:\Program Files\MATLAB71\work\dev\';
        num_samp = 32000;
        
        s(:,1) = wavread([Path 'female_src_2.wav'], num_samp); 
        s(:,2) = wavread([Path 'male_src_2.wav'], num_samp);
        s = s.';
        
        h = rand(3, 2);
        x = h * s;
        
        [z v] = pre_whitening(x);
        
        % Elimina os zeros
        minimo = min(min(abs(z(z ~= 0))));
        z(z == 0) = minimo;
        
        y = natICA(z(1:2, :), 'ScoreFunction', 'sign');
        
        figure; scatter3(z(1,:), z(2,:), z(3,:), 2, 'k');
        xlabel('z_1'); ylabel('z_2'); zlabel('z_3');
        
        figure; scatter(z(1,:), z(2,:), 2, 'k');
        set(gca, 'XTick', [], 'YTick', [], 'color', 'white'); xlabel('z_1'); ylabel('z_2');
        
        figure; scatter(z(1,:), z(3,:), 2, 'k');
        set(gca, 'XTick', [], 'YTick', [], 'color', 'white'); xlabel('z_1'); ylabel('z_3');
        
        figure; scatter(z(2,:), z(3,:), 2, 'k');
        set(gca, 'XTick', [], 'YTick', [], 'color', 'white'); xlabel('z_2'); ylabel('z_3');
        
end
