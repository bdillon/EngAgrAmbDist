#####################################################################
###                                                               ###
###             Agreement distance analysis pipeline              ###                                    ###
###                                                               ###
###                                                               ###
###             based on code by                                  ###
###             Christopher Hammerly                              ###
###             UMass Amherst - 04.22.16                          ###
###                                                               ###
#####################################################################

### Condition key for convenient reference: order is N1 N2 V; S = singular, P = plural
# agr-dist-a = grammatical, no intervener
# agr-dist-b = ungrammatical, no intervener
# agr-dist-c = grammatical, appositive intervener
# agr-dist-d = ungrammatical, appositive intervener
# agr-dist-e = grammatical, restrictive intervener
# agr-dist-f = ungrammatical, restrictive intervener

#### Load libraries

library(MASS)
library(dplyr)
library(ggplot2)
library(downloader)
library(ez)
library(gridExtra)

#### Set working directory
#### Eventually, I will add a github link so people will be able to download the data directly and run analyses seamlessly.
#### Here is a template for doing that:
#url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
#filename <- "msleep_ggplot2.csv"
#if (!file.exists(filename)) download(url,filename)

setwd("~/Dropbox/Current Work/EngAgrAmbDist/Figures")

#### Load data

mycols <- c("Subject","MD5","TrialType","Number","Element","Experiment","Item", "question", "response","null","RT","null")
results <- read.csv('~/Dropbox/Current Work/EngAgrAmbDist/Data/results.csv',
                    header = 0, 
                    sep = ",",
                    quote = "",
                    comment.char = "#",
                    col.names=mycols,
                    fill = TRUE)

#### segregate agreement distance data from fillers and other experiment

data <- subset(results , Experiment %in% c('agr-dist-a','agr-dist-b','agr-dist-c','agr-dist-d','agr-dist-e','agr-dist-f'))

### Adding columns for the factors: gram (+gram, -gram) by intervener type (none, appos, restr)

data$gram <- ifelse(data$Experiment %in% c('agr-dist-a','agr-dist-c','agr-dist-e'),'+GRAM','-GRAM')
data$intervener <- ifelse(
  data$Experiment %in% c('agr-dist-a','agr-dist-b'),
    'NONE',
    ifelse(data$Experiment %in% c('agr-dist-c','agr-dist-d'),'APPOS','RESTR')
  )

### Fixing preferred order of factors for purposes of graphing later

data$intervener <- factor(data$intervener, levels = c("NONE","RESTR","APPOS"))
data$gram <- factor(data$gram, levels = c("+GRAM","-GRAM"))

#### segregate acceptability judgments and reading time data

data.acceptability <- subset(data, TrialType == 'Question')

data.RT <- subset(data, TrialType == 'DashedSentence')
names(data.RT) <- c("Subject","MD5","TrialType","Number","Element","Experiment","Item", "region", "fragment","RT","null","sentence", "gram", "intervener")

#######################################
###                                 ###
###     Acceptability analysis      ###
###                                 ###
#######################################

#### a sanity check to ensure the correct number of subjects have been processed through
print(data.acceptability %>% summarise(number = n_distinct(Subject)))

#### change the values for the judgments from characters to numbers 

data.acceptability$response <- as.numeric(as.character(data.acceptability$response))
data.acceptability$z <- ave(data.acceptability$response, data.acceptability$Subject, FUN = scale)

#### a sanity check on distribution of judgments across experiment
pdf('agr-dist-alljudgments.pdf')
ggplot(data.acceptability,aes(x=response))+geom_histogram(binwidth=1)+xlim(1,7)+ggtitle("JUDGMENT DISTRIBUTION")
dev.off()


#### Basic descriptive stats
#### By condition means / SEM, for raw and scaled judgments

subj.by.cond <- data.acceptability %>% 
    group_by(Subject, gram,intervener) %>% 
    summarise(mean.resp = mean(response), 
              mean.z = mean(z))

cond.summ <- subj.by.cond %>%
  group_by(gram,intervener) %>%
  summarise(mean_cond = mean(mean.resp),
            SEM = sd(mean.resp)/sqrt(n_distinct(Subject)),
            average.z = mean(mean.z))

#### Basic descriptive stats
#### By intervener difference scores, for raw and scaled judgments, incl. 95%CI

differences.by.subj <- data.acceptability %>% 
    group_by(Subject,intervener) %>%
    filter(gram=='+GRAM') %>%
    summarise(gram.mean = mean(response),
              gram.z.mean = mean(z)
    )

differences.by.subj <- data.acceptability %>%
  group_by(Subject,intervener) %>%
  filter(gram=='-GRAM') %>%
    summarise(ungram.mean = mean(response),
              ungram.z.mean = mean(z)
    ) %>%
  right_join(differences.by.subj)

differences.summary <- differences.by.subj                                  %>%
  group_by(intervener)                                                      %>%
  summarise(
            N = n_distinct(Subject),
            mean.diff = mean(gram.mean-ungram.mean),
            mean.sem = sd(gram.mean-ungram.mean)/sqrt(n_distinct(Subject)),
            z.diff = mean(gram.z.mean-ungram.z.mean),
            z.sem = sd(gram.z.mean-ungram.z.mean)/sqrt(n_distinct(Subject))) %>%
  mutate(ci.lower.mean = mean.diff - qt(.975,df=N-1)*mean.sem,
         ci.upper.mean = mean.diff + qt(.975,df=N-1)*mean.sem,
         ci.lower.z = z.diff - qt(.975,df=N-1)*z.sem,
         ci.upper.z = z.diff + qt(.975,df=N-1)*z.sem
    )

#### Visualize data

 diff.plot <- ggplot(differences.summary,aes(x=intervener,y=mean.diff,fill=intervener))+geom_bar(position='dodge',stat = "identity")+geom_linerange(stat='identity',mapping=aes(ymax = ci.lower.mean,ymin=ci.upper.mean,width=.25))+ scale_fill_manual(values=c("navy","maroon","grey"))+xlab("INTERVENER TYPE")+ylab("GRAMMATICALITY PENALTY")
 raw.plot <- ggplot(cond.summ,aes(x=intervener,y=mean_cond,fill=gram))+geom_bar(position='dodge',stat = "identity")+geom_linerange(stat='identity',position = position_dodge(width = 0.9),mapping=aes(ymax = mean_cond+SEM,ymin=mean_cond-SEM)) + scale_fill_manual(values=c("navy","maroon"))+xlab("INTERVENER TYPE")+ylab("MEAN RATING")+ylim(0,7)

pdf('agr-dist-judgmentplot.pdf')
grid.arrange(raw.plot, diff.plot,ncol=2)
dev.off()

#### Hypothesis testing

ezANOVA(data = data.acceptability, dv = response, wid = Subject, within = c(gram, intervener))

###    
###    $ANOVA
###               Effect DFn DFd         F            p p<.05        ges
###    2            gram   1  47 30.271840 1.513651e-06     * 0.14515451
###    3      intervener   2  94  6.603727 2.072216e-03     * 0.01531282
###    4 gram:intervener   2  94  5.328157 6.427377e-03     * 0.01344144       
###    Sig interaction

###    $`Mauchly's Test for Sphericity`
###               Effect         W         p p<.05
###    3      intervener 0.9617592 0.4078737      
###    4 gram:intervener 0.9695211 0.4907012      ###    

###    $`Sphericity Corrections`
###               Effect       GGe       p[GG] p[GG]<.05      HFe       p[HF] p[HF]<.05
###    3      intervener 0.9631677 0.002375337         * 1.003514 0.002072216         *
###    4 gram:intervener 0.9704226 0.006969256         * 1.011565 0.006427377         *###    
###    Survives sphericity correction

#######################################
###                                 ###
###     Reading time analysis       ###
###                                 ###
#######################################

#### Sanity check for subject number

print(data.RT %>% summarise(number = n_distinct(Subject)))

#### Ensure RT data is numeric

data.RT$RT <- as.numeric(as.character(data.RT$RT))

#### Remove overlong outliers

data.RT <- filter(data.RT, RT < 5000)

#### Identify best transform for data; log tranform is winner

boxcox(RT~Experiment,data=data.RT)

#### Descriptive stats

RT.subj.by.cond <- data.RT %>%
  group_by(Subject, gram, intervener, region) %>%
  summarise(raw.rt = mean(RT),
            log.rt = mean(log(RT)))

RT.cond.summ <- RT.subj.by.cond %>%
  group_by(gram, intervener, region) %>%
  summarise(mean.raw = mean(raw.rt),
            raw.se = sd(raw.rt)/sqrt(n_distinct(Subject)),
            mean.log = mean(log.rt),
            log.se = sd(log.rt)/sqrt(n_distinct(Subject)))

#### 

ezANOVA(data = filter(data.RT, region == 4), dv = log(RT), wid = Subject, within = c(gram, intervener))

###### $ANOVA
######            Effect DFn DFd          F            p p<.05          ges
###### 2            gram   1  47 24.7337377 9.206386e-06     * 0.0199605363
###### 3      intervener   2  94  1.5205211 2.239252e-01       0.0018296481
###### 4 gram:intervener   2  94  0.1946658 8.234407e-01       0.0002722985###### 

###### $`Mauchly's Test for Sphericity`
######            Effect         W          p p<.05
###### 3      intervener 0.8218206 0.01096153     *
###### 4 gram:intervener 0.9235877 0.16069308      ###### 

###### $`Sphericity Corrections`
######            Effect       GGe     p[GG] p[GG]<.05       HFe     p[HF] p[HF]<.05
###### 3      intervener 0.8487671 0.2262212           0.8772331 0.2258793          
###### 4 gram:intervener 0.9290120 0.8076579           0.9656772 0.8160160          

#### 

ezANOVA(data = filter(data.RT, region == 4 & intervener != 'NONE'), dv = log(RT), wid = Subject, within = c(gram, intervener))

###### $ANOVA
######            Effect DFn DFd          F            p p<.05          ges
###### 2            gram   1  47 24.7337377 9.206386e-06     * 0.0199605363
###### 3      intervener   2  94  1.5205211 2.239252e-01       0.0018296481
###### 4 gram:intervener   2  94  0.1946658 8.234407e-01       0.0002722985###### 

###### $`Mauchly's Test for Sphericity`
######            Effect         W          p p<.05
###### 3      intervener 0.8218206 0.01096153     *
###### 4 gram:intervener 0.9235877 0.16069308      ###### 

###### $`Sphericity Corrections`
######            Effect       GGe     p[GG] p[GG]<.05       HFe     p[HF] p[HF]<.05
###### 3      intervener 0.8487671 0.2262212           0.8772331 0.2258793          
###### 4 gram:intervener 0.9290120 0.8076579           0.9656772 0.8160160          

#### Plotting RT data

pdf('agr-dist-RTplot.pdf',width=8.5,height=5)

#### Doing this for multiple regions; using for loop to avoid repetition

region.plots <- {}
for (cur.region in c('2','3','4')){
  region.plots[[cur.region]] <- ggplot(filter(RT.cond.summ,region==cur.region),aes(x=intervener,y=mean.raw,fill=gram))+geom_bar(position='dodge',stat = "identity")+geom_linerange(stat='identity',position = position_dodge(width = 0.9),mapping=aes(ymax = mean.raw+raw.se, ymin=mean.raw-raw.se)) + scale_fill_manual(values=c("navy","maroon"))+xlab("INTERVENER TYPE") + ylab("MEAN RT") + ggtitle(paste("RTs @ REGION ",cur.region,sep=""))+guides(fill=FALSE)
}
grid.arrange(region.plots[['2']],region.plots[['3']],region.plots[['4']],ncol=3)
dev.off()