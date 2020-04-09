theta = [10 10]; lob = [1e-1 1e-1]; upb = [20 20];
[dmodel, perf] = dacefit(S, Y, @regpoly0, @corrgauss, theta, lob, upb);
