# NFL-Play-Call-Predictions
R-Code and data for predicting play-calls in the National Football League (NFL)

Using the R-Code files, hidden Markov models are fitted to NFL play-by-play data (seasons 2009 - 2016) and 
play-calls of the 2017 season are predicted. The file "FullModelCode.R" fits a model to the whole data, i.e. one model
for all teams. The file "TeamModelCode.R" allows to fit a model for a chosen team. For both files, both the training and
test data are included in the .RData files "TrainingDataList.RData" and "TestDataList.RData", respectively.
