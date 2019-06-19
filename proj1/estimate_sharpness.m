function sharpness = estimate_sharpness(G)
    [Gx, Gy] = gradient(G);
    S = sqrt(Gx.*Gx + Gy.*Gy);
    sharpness=sum(sum(S))./(numel(Gx));
end