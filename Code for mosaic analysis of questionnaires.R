# Gender Mosaic Code for Questionnaires

# Read me ----
# This code excepts data files of the following structure: 
#The first column should be named "Group" and specify "F" for female, and "M" for male participants. 
#Besides the first column, you can include your questionnaire answers data that will be included in the mosaic calculation. Notice that this analysis will make sense only for questionnaires with ordinal answers.

#The following 3 functions will help you conduct the gender mosaic analysis: 

#1.Cohens_D_dist(Data):
#Input variable: "Data" - your data (should follow the instructions above).
#Result: this function returns a view of the Cohen's D values distribution and a dataframe with the Cohen's D values.


#2.General_Mosaic(Data,CohensD_Min,Percentage=0.33): generates the gender mosaic figure from a given data set.
#Input variables:
#"Data": your data (should follow the instructions above).
#"CohensD_Min": the minimum cohen’s D absolute value (has to be at least 0.5) of the variables that will be included in the mosaic.
#"Percentage" (Default = 0.33): The percentage that would be used to determine the "Female-end" (F) and "Male-end" (M) range of scores for each variable.

#Results:
# On the basis of the minimum cohen’s D value and percentage, this function returns one mosaic figure for males and one for females.
# Stats table, contains the following data: Variables name,
#Cohens_D,
#The response option which serves the limit of the Male-end,
#Percent of males with M,
# The response option which serves the limit of the Female-end,
#Percent of females with F,
#Percent of males with F,
#Percent of females with M 
# Notice: there are 4 colors in each mosaic : Blue for Male, Pink for Female, White for intermediate and Grey for missing values.

#3.Internally_consistent_function(Mosaic_M,Mosaic_F,allow_nas=0):  
#Input variables: a."Mosaic_M" and "Mosaic_F" that was returned from the General_Mosaic function (the mosaic tables that were used to create the figure).
#b."allow_nas"(Default = 0): number of missing values allowed for a row (observation) to still be considered internally consistent.
# Result: a.number of internally consistent males and females, meaning participants that all of their variables values are the same – all "M" or all "F"
#b.number of mosaic males and females, meaning participants that has at least one of their variables values classified as "M" and at least one classified as "F"


#Find Male-end and Female-end regions algorithm:
#1.If males’ average is larger than females’ average:
#Start from the highest answer towards the lowest. For each answer add the percentage of males that choose this answer to the cumulative sum, until you first reach the 
#input percentage. Calculate the difference of the cumulative sum before and after we added the last answer, from the input percentage.
#If the difference between the two differences is larger than 5(%), chose the option (include or not the last answer) that gave the smaller difference from the input percentage as the limit of the Male-end. 
#Otherwise, wait until we finish with the females’ region and then choose the option that is closest to the females’ actual percentage (if we have the same situation in the female side, we will compare all 4 options – 2 for females and 2 for males – and chose the one with the closest percentages).
#For the Female-end region, repeat the same procedure, only now we’ll start from the lowest answer number towards the highest.
#2. If males’ average is smaller than females’ average: 
#We’ll follow the same steps as in (1), only now we’ll start from the lowest answer number for the males’ region and from the highest answer number for the females.


#libraries----
library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)
library(dplyr)
library(lattice)
library(lsr)

#Functions----
Cohens_D_dist = function(Data){
  # Setup
  Number_of_Variables = ncol(Data) - 1
  Variables_Names = colnames(Data[,-1]) 
  
  # Creation of subgroups 
  Data_total = Data[,-1] #Data without Group  
  Data_M = Data_total[Data[,'Group']=="M",] #Males data
  Data_F = Data_total[Data[,'Group']=="F",] #Females data
  
  #Convert all data to numeric values
  Data_total <- apply(Data_total, 2, as.numeric)
  Data_M <- apply(Data_M, 2, as.numeric)
  Data_F <- apply(Data_F, 2, as.numeric)
  
  #calculating Cohen's d value for each variable
  Cohens_D=rep(NA,length(Number_of_Variables))
  for (i in 1:Number_of_Variables){
    sign = sign(mean(Data_F[,i], na.rm = T) - mean(Data_M[,i], na.rm = T))
    Cohens_D[i]=round(cohensD(Data_M[,i],Data_F[,i])*sign,3)
  }
  #plotting the Cohen's d distribution
  Cohen_D_df = data.frame(Variables_Name=Variables_Names ,  Cohens_D_value=Cohens_D)
  plot = ggplot(Cohen_D_df,aes(x=Cohens_D_value)) + geom_histogram(binwidth=.2,color="black", fill="white")+
    ggtitle("Cohen's d Distribution") + xlab("Cohen's d") + ylab("Count")
  
  
  return(list(plot,Cohen_D_df))
}

#----

General_Mosaic = function(Data,CohensD_Min,Percentage=0.33){
  
  # Setup
  Number_of_Variables = ncol(Data) - 1
  Variables_Names = colnames(Data[,-1]) 
  
  # Creation of subgroups 
  Data_total = Data[,-1] #Data without Group  
  Data_M = Data_total[Data[,'Group']=="M",] #Males data
  Data_F = Data_total[Data[,'Group']=="F",] #Females data
  
  #Convert all data to numeric values
  Data_total <- apply(Data_total, 2, as.numeric)
  Data_M <- apply(Data_M, 2, as.numeric)
  Data_F <- apply(Data_F, 2, as.numeric)
  
  # Calculate absolute value of the Cohen's D for the differences between Males and Females 
  Cohens_D=rep(NA,length(Number_of_Variables))
  cd_signs = rep(NA,length(Number_of_Variables))
  for (i in 1:Number_of_Variables){
    cd_signs[i] = sign(mean(Data_F[,i], na.rm = T) - mean(Data_M[,i], na.rm = T))
    Cohens_D[i]=round(abs(cohensD(Data_M[,i],Data_F[,i])),3)}
  
  #check the Cohen's D input value
  if (CohensD_Min<0.5){
    return("Selected Cohen's D is smaller than 0.5, please select a bigger one")}
  
  if (max(Cohens_D)<CohensD_Min){
    return("Selected Cohen's D is too big, please select a smaller one")}
  
  #filtering variables with small Cohen's D
  relevant_var = c(which(Cohens_D>CohensD_Min))
  Data_total = Data_total[,relevant_var]
  Data_M = Data_M[,relevant_var]
  Data_F = Data_F[,relevant_var]
  Number_of_Variables = length(relevant_var)
  Variables_Names = Variables_Names[relevant_var]
  
  
  # Calculation of average values of females and Males for each variable
  males_mean = apply(Data_M,2,mean,na.rm = TRUE)
  females_mean = apply(Data_F,2,mean,na.rm = TRUE)
  
  
  #Calculate the values for the males' and females' region limits
  males_region_limits = vector(mode = 'numeric', length = Number_of_Variables)
  females_region_limits = vector(mode = 'numeric', length = Number_of_Variables)
  males_region_limits_per = vector(mode = 'numeric', length = Number_of_Variables)
  females_region_limits_per = vector(mode = 'numeric', length = Number_of_Variables)
  
  for (i in c(1:Number_of_Variables)) {
    males_Percentage_sum <- 0 
    females_Percentage_sum <- 0 
    variable_data_without_na <- na.omit(as.data.frame(list(Data_total[,i], Data[,1])))
    prc_males <- round((table(variable_data_without_na[,1],variable_data_without_na[,2])[,'M'] / sum(variable_data_without_na[,2]=='M'))*100 ,2)
    prc_females <- round((table(variable_data_without_na[,1],variable_data_without_na[,2])[,'F'] / sum(variable_data_without_na[,2]=='F'))*100 ,2)
    
 #----------------   
    if(males_mean[i]>females_mean[i]){ 
    #finding males limit:
      wait_for_females_limit <- FALSE
      males_limit <- max(variable_data_without_na[,1])
      if (prc_males[as.character(males_limit)] >= Percentage*100){ # if the first answer exceeds the given percent
        males_region_limits[i] <- males_limit
        males_region_limits_per[i] <- prc_males[as.character(males_limit)]
      }

      else{
        while (males_Percentage_sum<Percentage*100) { 
          males_Percentage_sum <- males_Percentage_sum + prc_males[as.character(males_limit)]
          males_limit <- males_limit-1}
        
        if ((males_Percentage_sum - (Percentage*100)) < ((Percentage*100) - (males_Percentage_sum - prc_males[as.character(males_limit+1)]))-5){
          males_region_limits[i] <- (males_limit+1)
          males_region_limits_per[i] <- males_Percentage_sum
        }
        else if ((males_Percentage_sum - (Percentage*100)) > ((Percentage*100) - (males_Percentage_sum - prc_males[as.character(males_limit+1)]))+5){
          males_region_limits[i] <- (males_limit+2)
          males_region_limits_per[i] <- males_Percentage_sum - prc_males[as.character(males_limit+1)]
        }
        else{
          wait_for_females_limit<-TRUE}
      }
      
    #finding females limit:
      wait_for_males_limit <- FALSE
      females_limit <- min(variable_data_without_na[,1])
      if (prc_females[as.character(females_limit)] >= Percentage*100){ # if the first answer exceeds the given percent
        females_region_limits[i] <- females_limit
        females_region_limits_per[i] <- prc_females[as.character(females_limit)]
      }
      else{
        while (females_Percentage_sum<Percentage*100) { 
          females_Percentage_sum <- females_Percentage_sum + prc_females[as.character(females_limit)]
          females_limit <- females_limit+1}
        
        if ((females_Percentage_sum - (Percentage*100)) < ((Percentage*100) - (females_Percentage_sum - prc_females[as.character(females_limit-1)]))-5){
          females_region_limits[i] <- (females_limit-1)
          females_region_limits_per[i] <- females_Percentage_sum
        }
        else if ((females_Percentage_sum - (Percentage*100)) > ((Percentage*100) - (females_Percentage_sum - prc_females[as.character(females_limit-1)]))+5){
          females_region_limits[i] <- (females_limit-2)
          females_region_limits_per[i] <- females_Percentage_sum - prc_females[as.character(females_limit-1)]
        }
        else{
          wait_for_males_limit<-TRUE
        }}  
      
      #breaking ties
      if (wait_for_females_limit&(!wait_for_males_limit)){
        if (abs(males_Percentage_sum - females_region_limits_per[i])<=
        abs((males_Percentage_sum - prc_males[as.character(males_limit+1)])-females_region_limits_per[i])){
          males_region_limits[i] <- (males_limit+1)
          males_region_limits_per[i] <- males_Percentage_sum}
        else{
          males_region_limits[i] <- (males_limit+2)
          males_region_limits_per[i] <- males_Percentage_sum - prc_males[as.character(males_limit+1)]
        }}
      else if (wait_for_males_limit&(!wait_for_females_limit)){
        if (abs(females_Percentage_sum - males_region_limits_per[i])<=
            abs((females_Percentage_sum - prc_females[as.character(females_limit-1)])-males_region_limits_per[i])){
          females_region_limits[i] <- (females_limit-1)
          females_region_limits_per[i] <- females_Percentage_sum}
          
        else{
          females_region_limits[i] <- (females_limit-2)
          females_region_limits_per[i] <- females_Percentage_sum - prc_females[as.character(females_limit-1)]}
      }
      
      else if (wait_for_males_limit&wait_for_females_limit){
        males_options <-list(c(males_Percentage_sum, males_Percentage_sum - prc_males[as.character(males_limit+1)]),
                             c(males_limit+1,males_limit+2)) 
        females_options <- list(c(females_Percentage_sum,females_Percentage_sum - prc_females[as.character(females_limit-1)]),
                               c(females_limit-1,females_limit-2))
        
        index <- which.min(
                  c(abs(males_options[[1]][1]-females_options[[1]][1]),
                  abs(males_options[[1]][1]-females_options[[1]][2]),
                  abs(males_options[[1]][2]-females_options[[1]][1]),
                  abs(males_options[[1]][2]-females_options[[1]][2])))
        
        options_code <- list(c(1,1),c(1,2),c(2,1),c(2,2))
  
        males_region_limits[i] <- males_options[[2]][options_code[[index]][1]]
        males_region_limits_per[i] <- males_options[[1]][options_code[[index]][1]]
        females_region_limits[i] <- females_options[[2]][options_code[[index]][2]]
        females_region_limits_per[i] <-females_options[[1]][options_code[[index]][2]]}
      if (males_region_limits[i]<=females_region_limits[i]){
        return("Warning! Regions limits intersect, try a smaller percentage")}
    }
    
    
#----------------   
    if(males_mean[i]<females_mean[i]){ 
      #finding males limit:
      wait_for_females_limit <- FALSE
      males_limit <- min(variable_data_without_na[,1])
      if (prc_males[as.character(males_limit)] >= Percentage*100){ # if the first answer exceeds the given percent
        males_region_limits[i] <- males_limit
        males_region_limits_per[i] <- prc_males[as.character(males_limit)]
      }
      
      else{
        while (males_Percentage_sum<Percentage*100) { 
          males_Percentage_sum <- males_Percentage_sum + prc_males[as.character(males_limit)]
          males_limit <- males_limit+1}
        
        if ((males_Percentage_sum - (Percentage*100)) < ((Percentage*100) - (males_Percentage_sum - prc_males[as.character(males_limit-1)]))-5){
          males_region_limits[i] <- (males_limit-1)
          males_region_limits_per[i] <- males_Percentage_sum
        }
        else if ((males_Percentage_sum - (Percentage*100)) > ((Percentage*100) - (males_Percentage_sum - prc_males[as.character(males_limit-1)]))+5){
          males_region_limits[i] <- (males_limit-2)
          males_region_limits_per[i] <- males_Percentage_sum - prc_males[as.character(males_limit-1)]
        }
        else{
          wait_for_females_limit<-TRUE}
      }
      
      #finding females limit:
      wait_for_males_limit <- FALSE
      females_limit <- max(variable_data_without_na[,1])
      if (prc_females[as.character(females_limit)] >= Percentage*100){ # if the first answer exceeds the given percent
        females_region_limits[i] <- females_limit
        females_region_limits_per[i] <- prc_females[as.character(females_limit)]
      }
      else{
        while (females_Percentage_sum<Percentage*100) { 
          females_Percentage_sum <- females_Percentage_sum + prc_females[as.character(females_limit)]
          females_limit <- females_limit-1}
        
        if ((females_Percentage_sum - (Percentage*100)) < ((Percentage*100) - (females_Percentage_sum - prc_females[as.character(females_limit+1)]))-5){
          females_region_limits[i] <- (females_limit+1)
          females_region_limits_per[i] <- females_Percentage_sum
        }
        else if ((females_Percentage_sum - (Percentage*100)) > ((Percentage*100) - (females_Percentage_sum - prc_females[as.character(females_limit+1)]))+5){
          females_region_limits[i] <- (females_limit+2)
          females_region_limits_per[i] <- females_Percentage_sum - prc_females[as.character(females_limit+1)]
        }
        else{
          wait_for_males_limit<-TRUE
        }}  
      
      #breaking ties
      if (wait_for_females_limit&(!wait_for_males_limit)){
        if (abs(males_Percentage_sum - females_region_limits_per[i])<=
            abs((males_Percentage_sum - prc_males[as.character(males_limit-1)])-females_region_limits_per[i])){
          males_region_limits[i] <- (males_limit-1)
          males_region_limits_per[i] <- males_Percentage_sum}
        else{
          males_region_limits[i] <- (males_limit-2)
          males_region_limits_per[i] <- males_Percentage_sum - prc_males[as.character(males_limit-1)]
        }}
      else if (wait_for_males_limit&(!wait_for_females_limit)){
        if (abs(females_Percentage_sum - males_region_limits_per[i])<=
            abs((females_Percentage_sum - prc_females[as.character(females_limit+1)])-males_region_limits_per[i])){
          females_region_limits[i] <- (females_limit+1)
          females_region_limits_per[i] <- females_Percentage_sum}
        
        else{
          females_region_limits[i] <- (females_limit+2)
          females_region_limits_per[i] <- females_Percentage_sum - prc_females[as.character(females_limit+1)]}
      }
      
      else if (wait_for_males_limit&wait_for_females_limit){
        males_options <-list(c(males_Percentage_sum, males_Percentage_sum - prc_males[as.character(males_limit-1)]),
                             c(males_limit-1,males_limit-2)) 
        females_options <- list(c(females_Percentage_sum,females_Percentage_sum - prc_females[as.character(females_limit+1)]),
                                c(females_limit+1,females_limit+2))
        
        index <- which.min(
          c(abs(males_options[[1]][1]-females_options[[1]][1]),
            abs(males_options[[1]][1]-females_options[[1]][2]),
            abs(males_options[[1]][2]-females_options[[1]][1]),
            abs(males_options[[1]][2]-females_options[[1]][2])))
        
        options_code <- list(c(1,1),c(1,2),c(2,1),c(2,2))
        
        males_region_limits[i] <- males_options[[2]][options_code[[index]][1]]
        males_region_limits_per[i] <- males_options[[1]][options_code[[index]][1]]
        females_region_limits[i] <- females_options[[2]][options_code[[index]][2]]
        females_region_limits_per[i] <-females_options[[1]][options_code[[index]][2]]}
      
      if (males_region_limits[i]>=females_region_limits[i]){
        return("Warning! Regions limits intersect, try a smaller percentage")}
    }}
  
#----------------    
  
  # Creating the Mosaic 
  percent_of_females_classified_as_males <- vector(mode = 'numeric', length = Number_of_Variables)
  percent_of_males_classified_as_females <- vector(mode = 'numeric', length = Number_of_Variables)  
  Mosaic <- as.data.frame(matrix(NA,nrow = nrow(Data_total),
                               ncol = (Number_of_Variables+1)))
  colnames(Mosaic)=c('Group',Variables_Names)
  Mosaic[,1]= Data$Group
  
  #Classifying the observations as -1=Female,1=Male,0=Neutral or missing data, according to the Male and Female limits
  for (i in c(1:Number_of_Variables)) {
    if(males_mean[i]>females_mean[i]){
      Mosaic[,i+1] = ifelse(Data_total[,i]>=males_region_limits[i],1,
                            ifelse(Data_total[,i]<=females_region_limits[i],-1,0))}
    if(females_mean[i]>males_mean[i]){
      Mosaic[,i+1] = ifelse(Data_total[,i]<=males_region_limits[i],1,
                            ifelse(Data_total[,i]>=females_region_limits[i],-1,0))}
    percent_of_females_classified_as_males[i] <- (sum(na.omit(Mosaic[Mosaic[,1]=='F',i+1])==1)/length(na.omit(Mosaic[Mosaic[,1]=='F',i+1])))*100
    percent_of_males_classified_as_females[i] <- (sum(na.omit(Mosaic[Mosaic[,1]=='M',i+1])==-1)/length(na.omit(Mosaic[Mosaic[,1]=='M',i+1])))*100
  }
  Mosaic[is.na(Mosaic)]<-(-2)
  Mosaic_stats <- cbind("Variables Name"=Variables_Names,
                        "Cohens_D"=Cohens_D[relevant_var]*cd_signs[relevant_var],
                        "Males region limit "=males_region_limits,
                        "Percent of males that cross males limit"=males_region_limits_per,
                        "Females region limit"=females_region_limits,
                        "Percent of females that cross females limit"=females_region_limits_per,
                        "Percent of males classified as females"=round(percent_of_males_classified_as_females,3),
                        "Percent of females classified as males"=round(percent_of_females_classified_as_males,3))
  
  #Preparing the data for the gender mosaic figure
  Mosaic_M = Mosaic[Mosaic$Group=="M",-1]
  Mosaic_F = Mosaic[Mosaic$Group=="F",-1]
  rownames(Mosaic_M) <- 1:nrow(Mosaic_M)
  rownames(Mosaic_F) <- 1:nrow(Mosaic_F)
  
  colors <- colorRampPalette(c("#CCCCCC","#FF99FF","white", "#0099FF"))
  #males mosaic
  M_M = levelplot(t(Mosaic_M[nrow(Mosaic_M):1,]),col.regions = colors,asp=2,at=seq(-2, 1, by = 0.75),cut=3,xlab="Variables",
                  main=list(label = "Male Mosaic", hjust = 0.8),
                  colorkey=list(height = 0.8, width=0.6, space="right",at=seq(-2, 1, by = 0.75),
                                labels=list(at=c(-1.6, -0.9, -0.1, 0.6),labels=c("NA","F", "Neutral",  "M"))),
                  scales=list(x=list(rot=60,cex=.5,alternating=2),y=list(alternating=2,draw=FALSE),tck = c(0,1)))
  #females mosaic
  M_F = levelplot(t(Mosaic_F[nrow(Mosaic_F):1,]),col.regions = colors,asp=2,at=seq(-2, 1, by = 0.75),
                  cut=3,xlab="Variables", main=list(label = "Female Mosaic", hjust = 0.8),
                  colorkey=list(height = 0.8, width=0.6, space="right",at=seq(-2, 1, by = 0.75),
                                labels=list(at=c(-1.6, -0.9, -0.1, 0.6),labels=c("NA","F", "Neutral", "M"))),
                  scales=list(x=list(rot=60,cex=.5,alternating=2),y=list(alternating=2,draw=FALSE),tck = c(0,1)))
  #combined plot  
  plot=grid.arrange(M_M, M_F, ncol=2, nrow = 1,top=textGrob(paste("Mosaic Using",  Percentage*100, "% as Cutoff", sep=" ")))
  return(list(plot,Mosaic_M,
              Mosaic_F,Mosaic_stats))} 

#---- 

Internally_consistent_function = function(Mosaic_M,Mosaic_F,allow_nas = 0){
  
  check_row_internally_consistent= function(row,allow_nas){ 
    #checks if a row is internally consistent, allowing "allow_nas" number of missing values  
    if (sum(row==-2)<=allow_nas){
      all(row[!(row==-2)] == 1) || all(row[!(row==-2)] == -1)}
    else {FALSE}
  }
  
  consistent_Female = sum(apply(Mosaic_F, 1, FUN=function(row) check_row_internally_consistent(row, allow_nas=allow_nas))) #all 'male' or all 'female'
  consistent_Male <- sum(apply(Mosaic_M, 1, FUN=function(row) check_row_internally_consistent(row, allow_nas=allow_nas))) #all 'male' or all 'female'
  Mosaic_check_Male <- sum(apply(Mosaic_M, 1, function(row) any(row == 1) && any(row == -1)),na.rm = TRUE) #some 'male' and some 'female'
  Mosaic_check_Female <- sum(apply(Mosaic_F, 1, function(row) any(row == 1) && any(row == -1)),na.rm = TRUE) #some 'male' and some 'female'
  
  return(c("Internally consistent Females:"=consistent_Female,
           " Internally consistent Males:"=consistent_Male,
           " Mosaic Female participants:"=Mosaic_check_Female,
           " Mosaic Male participants:"=Mosaic_check_Male))
}

#----

#Using the functions above the create your gender mosaic:

df <- #Your data
Cohens_D_dist(df) #Looking at the Cohen's D values

results <- General_Mosaic(Data=df,CohensD_Min=0.5,Percentage=0.33) #creating the mosaic
results[[1]] #if no figure was created,check if the function returned a warning
Mosaic_M = results[[2]]
Mosaic_F = results[[3]]
Mosaic_stats = as.data.frame(results[[4]]) #Get mosaic stats table

#mosaic summary: number of internally consistent and mosaic participants
Internally_consistent_function(Mosaic_M,Mosaic_F,allow_nas = 0)  







