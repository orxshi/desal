H = 0.1;
d = 0.1;
mx = 21;


pLL = newp; Point(pLL) = {-H,-H,-H}; // lower left
pUL = newp; Point(pUL) = {-H,H,-H};  // upper left
pLR = newp; Point(pLR) = {H,-H,-H};  // lower right
pUR = newp; Point(pUR) = {H,H,-H};   // upper right

lL = newl; Line(lL) = {pLL,pUL};
lB = newl; Line(lB) = {pLR,pLL};
lR = newl; Line(lR) = {pUR,pLR};
lT = newl; Line(lT) = {pUL,pUR};

Transfinite Line{lL} = mx;
Transfinite Line{lR} = mx;
Transfinite Line{lB} = mx;
Transfinite Line{lT} = mx;

Line Loop(1) = {lL, lT, lR, lB};
sBase = news; Plane Surface(sBase) = {1};

Transfinite Surface{sBase};
Recombine Surface{sBase};

out[] = Extrude {0,0,d} {Surface{sBase}; Layers{1}; Recombine;};

Physical Surface(2) = {out[2], out[4], out[5]};
Physical Surface(4) = {sBase, out[0]};
Physical Surface(3) = {out[3]};
Physical Volume(1) = {out[1]};
