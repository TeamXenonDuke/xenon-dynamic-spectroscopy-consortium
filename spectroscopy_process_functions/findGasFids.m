function gasFrames = findGasFids(fids)
    MAX = max(abs(fids));
    MAX2 = [MAX(2:end),MAX(end)];
    diff = abs(MAX-MAX2);
    [foo,I] = max(diff);
    gasFrames = size(fids,2) - I;
end
    