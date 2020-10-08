function decibels = PSNR(A,B)

% PSNR  Find the PSNR (peak signal-to-noise ratio) between two intensity images A and B, each having values in the interval  [0,1]. 
% The answer is in decibels (dB).

if A == B
   error('Images are identical: PSNR has infinite value')
end
max2_A = max(max(A));
max2_B = max(max(B));
min2_A = min(min(A));
min2_B = min(min(B));
if max2_A > 1 || max2_B > 1 || min2_A < 0 || min2_B < 0
   error('input matrixes must have values in the interval [0,1]')
end
error_diff = A - B;
decibels = 20*log10(1/(sqrt(mean(mean(error_diff.^2)))));
% disp(sprintf('PSNR = +%5.2f dB',decibels))

% if A == B
%    error('Images are identical: PSNR has infinite value')
% end
% max2_A = max(max(A));
% max2_B = max(max(B));
% min2_A = min(min(A));
% min2_B = min(min(B));
%
% if max2_A > 255 || max2_B > 255 || min2_A < 0 || min2_B < 0
%   error('input matrixes must have values in the interval [0,255]')
% end
% error_diff = A - B;
% decibels = 20*log10(255/(sqrt(mean(mean(error_diff.^2)))));
% disp(sprintf('PSNR = +%5.2f dB',decibels))

end
