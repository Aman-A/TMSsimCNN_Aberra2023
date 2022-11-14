function outdata = loglinlogtransform(indata, linlimit)

% linlimit: e.g., 0.01. Data within +/- linlimit will mapped to [-1, 1].
% Data outside of linlimit are log-transformed and shifted by +/-
% (log10(linlimit)+1).
% Written by Boshuo Wang
% Duke University, 2022
outdata = indata;

% scale = 2*(log10(loglimit)-log10(linlimit)) + 2 * linlimit;

ind_pos = (indata > linlimit);
outdata(ind_pos) = 	(log10( indata(ind_pos))-log10(linlimit))+1;

ind_neg = (indata < -linlimit);
outdata(ind_neg) = -(log10(-indata(ind_neg))-log10(linlimit))-1;

ind_lin = abs(indata)<=linlimit;
outdata(ind_lin) = indata(ind_lin)/linlimit;



end
