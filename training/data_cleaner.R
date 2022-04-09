
data_path <- 'C:/Users/parfait/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/Data/mergedUPdata_2018-10-22_12_12.csv'
data_path <- 'C:/Users/eparfmu/Box Sync/PhD file/projects/Survivor analysis/Paper 3 Software quality control/Data/mergedUPdata_2018-10-22_12_12.csv'

# Import the Data
UPdata<- read.csv(data_path) # Follow the promped window to find the dataset csv file  is stored and double click on it
TRdata<-data.frame(Date=as.Date(UPdata$Created),TRcount=UPdata$Internal.real.TRs+
                     UPdata$Real.External.TRs,ReleaseType=UPdata$Release.Type,
                   Commits=UPdata$Number.of.Commits,
                   ChangedModule=UPdata$Number.of.changed.modules,
                   CorrectedTR=UPdata$Final.File.correct.TR,
                   RowAdded=UPdata$Final.Rows.Added,
                   RowRemoved=UPdata$Final.Rows.Removed,
                   ncloc=UPdata$Final.ncloc,
                   racoamRepo=UPdata$Final.in.Repo.racoam,
                   Cfiles=UPdata$Final.Num.files.c..,
                   JavaFiles=UPdata$Final.Num.files.java,
                   AMfiles=UPdata$Final.Number.of.M.A.Files,
                   ScannedFiles=UPdata$Final.Files.We.Are.Scanning,
                   FileComplexity=UPdata$Final.File.Complexity)

# Select only track 6
TRdata<-TRdata[grep("CXP9024418/6_",UPdata$Upgrade.Package),]
TRdata<-TRdata[order(TRdata$Date),]
TRdata<-as.data.frame(TRdata)
# define the interval time
GA.Dates<-TRdata$Date[which(TRdata$ReleaseType=="GA")]

time.intervals<-seq(min(TRdata$Date),max(TRdata$Date),
                    length.out =as.numeric(floor(max(TRdata$Date)-min(TRdata$Date))/30)) # set the time intervals to month periods

time.intervals[1]<-min(TRdata$Date)-1
# clean the data
TRdata<-na.omit(TRdata) # remove NA
TRdata<-TRdata[which(TRdata$AMfiles>0),] # remove data points with zero number of files
#TRdata<-TRdata[which(TRdata$ScannedFiles>0),]
TRdata<-TRdata[which(TRdata$ncloc>0),]
# modify variables
TRdata$Commits<-TRdata$Commits/TRdata$AMfiles # commits per file
TRdata$CorrectedTR<-TRdata$CorrectedTR/TRdata$ncloc
TRdata$Cfiles<-TRdata$Cfiles/TRdata$AMfiles
TRdata$JavaFiles<-TRdata$JavaFiles/TRdata$AMfiles
TRdata<-TRdata[,c(1,2,4,5,6,10,11,12,15)]
# standardise the data
TRdata[,-c(1,2,6)]<-apply(TRdata[,-c(1,2,6)],2,function(x)log(x+1)) 
#TRdata[,-c(1,2,6)]<-apply(TRdata[,-c(1,2,6)],2,function(x)(x-min(x))/max(x))

windows(24,24)
pairs(TRdata[,c(-1,-6)],labels = c("TR","NC","CM","NFC","CF","JF","FC"),lower.panel = NULL,cex=.7)

#Test data and training data
Training.Last.Day<-as.Date("2018-07-27") # set the last day of training the model
Test.Last.Day<-(Training.Last.Day+1)+30 # set the last day of test data
# Test data
Test.data<-TRdata[which(TRdata$Date<=Test.Last.Day&TRdata$Date>Training.Last.Day),] # test data
# Training data
Data<-TRdata[which(TRdata$Date<=Training.Last.Day),] # train data
time.intervals<-time.intervals[-which(time.intervals>Training.Last.Day)]
if(max(time.intervals)!=Training.Last.Day){time.intervals<-c(time.intervals,Training.Last.Day)}


