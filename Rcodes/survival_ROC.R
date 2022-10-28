
library(survivalROC)
survivalROC(Stime=merge.train.pred$time, status=merge.train.pred$status,marker=merge.train.pred$Risk_Score,predict.time=5475,method="KM")$AUC
Mayo_Grade = survivalROC(Stime=merge.train.pred$time, status=merge.train.pred$status,marker=merge.train.pred$Tumor_Grade,predict.time=5475,method="KM")
Mayo_Age = survivalROC(Stime=merge.train.pred$time, status=merge.train.pred$status,marker=merge.train.pred$Age,predict.time=5475,method="KM")
Mayo_Gender = survivalROC(Stime=merge.train.pred$time, status=merge.train.pred$status,marker=merge.train.pred$Gender,predict.time=5475,method="KM")
Mayo_IDH_Mut_Status = survivalROC(Stime=merge.train.pred$time, status=merge.train.pred$status,marker=merge.train.pred$IDH_Mut_Status,predict.time=5475,method="KM")
Mayo_PQ = survivalROC(Stime=merge.train.pred$time, status=merge.train.pred$status,marker=merge.train.pred$PQ,predict.time=5475,method="KM")
Mayo

plot(Mayo$FP,Mayo$TP,type="l",xlim=c(0,1),ylim=c(0,1),xlab="FP",ylab="TP",col=c("#a50f15"),lwd=2)
lines(Mayo_IDH_Mut_Status$FP,Mayo_IDH_Mut_Status$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#de2d26"),lwd=2)
lines(Mayo_Age$FP,Mayo_Age$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#fb6a4a"),lwd=2)
lines(Mayo_Grade$FP,Mayo_Grade$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#fcae91"),lwd=2)
lines(Mayo_Gender$FP,Mayo_Gender$TP,type="l",xlim=c(0,1),ylim=c(0,1),col=c("#fee5d9"),lwd=2)

legend('bottomright', legend=c("Risk Score", "IDH Mut Status","Age","Grade","Gender"),col=c("#a50f15", "#de2d26","#fcae91","#fcae91","#fee5d9"), cex=0.5,lty=1,lwd=2)
abline(0,1,lwd=2)