function [SDR, SIR, SAR]=bss_eval_freq(se, s)
% BSS_EVAL_FREQ Ordering and measurement of the separation quality for
% estimated source signals in terms of filtered true source, interference
% and artifacts, on each frequency
%
% [SDR, SIR, SAR]=bss_eval_freq(se, s)
%
% Inputs:
% se: nfreqs x nsrc x nsampl matrix containing the frequency representation
%     of estimated sources
% s: nfreqs x nsrc x nsampl matrix containing the frequency representation
%    of true sources
%
% Outputs:
% SDR: nsrc x nfreqs vector of Signal to Distortion Ratios
% SIR: nsrc x nfreqs vector of Source to Interference Ratios
% SAR: nsrc x nfreqs vector of Sources to Artifacts Ratios
%
% Observation: The modulus of the frequency representation is taken !!!
%

K = size(se, 1);
N = size(se, 2);

SIR = zeros(N, K);
SDR = zeros(N, K);
SAR = zeros(N, K);

for ck = 1:K
    tmpse = abs(squeeze(se(ck, :, :)));
    tmps = abs(squeeze(s(ck, :, :)));
    [SDR(:, ck), SIR(:, ck), SAR(:, ck), dummy]=bss_eval_sources(tmpse, tmps);
end
