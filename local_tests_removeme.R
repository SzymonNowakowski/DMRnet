data(miete)
Xtr <- miete[1:1500,-1]
ytr <- miete[1:1500,1]
Xte <- miete[1501:2053,-1]

m1 <- DMR(Xtr, ytr)
