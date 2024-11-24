
function nn = norm_lr(X,Y)
    [QX,RX]=qr(X);
    [QY,RY]=qr(Y);
    nn = norm(RX*RY','fro');
end


%function U = HestonExplicitClassicCNXYRC01(params,K,r,q,S,V,T,mode)