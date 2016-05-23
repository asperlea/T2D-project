library(randomForest)

print("Reading in clinical data . . .")
clinData=na.omit(read.table("Clinical data.csv", head=TRUE, sep=','))
print ("Done.")
binaryStatus=clinData[which(clinData$F_GLTOL5_6 != '?'),]
binaryStatus$F_GLTOL5_6=as.numeric(binaryStatus$F_GLTOL5_6)
binaryStatus[which(binaryStatus$F_GLTOL5_6 != 5), "F_GLTOL5_6"] <- 0
binaryStatus[which(binaryStatus$F_GLTOL5_6 == 5), "F_GLTOL5_6"] <- 1

print("Reading in methilation data . . .")
metData=read.table("metsimSample.txt", head=TRUE, row.names = 1, sep='\t')
print("Done.")
IDs=c('2273', '2320', '2368', '2403', '2421', '2436', '2469', '2480', '2481', '1796', '1813', '1944', '2129', '2226', '2230', '2262', '2264', '1208', '1234', '1280', '1301', '1436', '1517', '1530', '1738', '1744', '718', '758', '787', '816', '1004', '1049', '1058', '1082', '1147', '377', '379', '441', '495', '6564', '6632', '6639', '6807', '535', '539', '634', '635', '697', '10886', '10280', '10377', '10390', '10437', '10406', '10548', '10567', '10664', '10879', '7366', '7432', '7604', '7626', '10012', '10016', '10039', '10070', '10071', '1148', '1452', '1545', '1674', '1768', '1896', '2057', '2348', '1559', '1165', '1678', '1455', '5521', '1777', '2100', '1932', '6854', '6875', '7096', '7136', '2430', '5530', '1264', '1457', '1580', '1702', '1779', '1942', '2103', '2529', '5571', '1303', '1462', '1583', '1722', '1822', '1974', '2176', '2655', '5812', '1409', '1483', '1590', '1725', '1840', '1984', '2242', '2656', '5833', '1438', '1490', '1599', '1736', '1841', '1987', '2261', '2760', '5849', '1445', '1502', '7295', '5242', '5283', '5365', '1619', '1747', '1885', '1999', '2267', '5273', '5859', '1446', '1518', '1653', '1757', '1891', '2290', '2044', '5462', '5919', '1449', '1524', '1666', '1764', '1892', '2048', '2321', '5470', '6194', '6514', '6515', '6537', '6648', '6666', '6684', '6761', '5922', '6083', '6101', '6183', '6200', '6203', '6277', '6283', '5392', '5481', '5668', '5949', '6425', '5923', '5935', '5951', '5953', '5954', '6005', '6010', '6028', '6041', '6106', '6501', '4479', '4481', '4509', '4518', '4544', '4856', '5084', '5166', '5210', '2513', '2726', '2741', '2743', '2748', '2808', '2902', '4008', '4219', '10177', '10419', '10896', '1111', '2680', '4383', '4410', '520', '522', '5263', '5394', '5460', '569', '5790', '582', '6311', '6737', '6978', '7742', '8219', '844', '845', '874')
#withheldIDs=c('787', '1264', '4219', '5242', '6311', '6978', '7742')
colnames(metData) <- IDs
metData<-t(metData)
metData<-metData[as.character(sort(as.numeric(rownames(metData)))),]

diseaseStatus=c()
for (ID in sort(as.numeric(IDs))) {
  if (ID %in% binaryStatus$ID) {
    diseaseStatus=c(diseaseStatus, binaryStatus$F_GLTOL5_6[which(binaryStatus$ID == ID)])
  } else {
    #print(which(binaryStatus$ID == ID))
    #print(length(cbind(metData[which(rownames(metData) == ID),], binaryStatus$F_GLTOL5_6[which(binaryStatus$ID == ID)])))
    metData<-metData[-which(rownames(metData) == ID),]
  }
}

metData<-cbind(metData, diseaseStatus)
metData<-data.frame(metData)
metData$diseaseStatus<-as.factor(metData$diseaseStatus)

print("Imputing missing values in data . . .")
metDataImputed<-rfImpute(diseaseStatus ~ ., metData)
print("Done.")

#metDataImputed<-metData
avgAccuracy<-0
for (i in 1:10) {
  print("Sampling testing data . . .")
  testingData<-metDataImputed[sample(nrow(metDataImputed), 5), ]
  print("Done.")
  metDataImputedCopy<-metDataImputed

  for (ID in rownames(testingData)) {
    pos<-which(rownames(metDataImputedCopy) == ID)
    metDataImputedCopy<-metDataImputedCopy[-pos,]
  }
  
  print("Fitting model . . .")
  fit <- randomForest(diseaseStatus ~ ., data=metDataImputedCopy)
  print("Done.")
  
  prediction <- predict(fit, testingData)
  TP<-0
  TN<-0
  FP<-0
  FN<-0
  for (i in 1:length(prediction)) {
    if (testingData$diseaseStatus[i] == 1) {
      if (prediction[i] == 1) {
        TP<-TP+1
      } else {
        FN<-FN+1
      }
    } else {
      if (prediction[i] == 0) {
        TN<-TN+1
      } else {
        FP<-FP+1
      }
    }
  }
  print(paste(TP, FP, TN, FN))
  avgAccuracy<-avgAccuracy+((TP+TN)/(TP+FP+FN+TN))
}
print(paste("Average accuracy", avgAccuracy / 10))
