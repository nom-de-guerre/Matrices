Data <- read.table ("data")
T <- Data[Data$V1 == "DJS_RESULT", ]
T <- T[order (T$V2), ]
E <- Data[Data$V1 == "DJS_INFER", ]
plot (T$V3 ~ T$V2, col = "blue", pch=19); lines (T$V3 ~ T$V2, col = "blue")
points (T$V4 ~ T$V2, col="red")
points (E$V4 ~ E$V2, pch=21, col="orange")

