#install packages
install.packages("tm")
install.packages("slam")
install.packages("topicmodels")

#load library
library("tm")
library("slam")
library("topicmodels")

#set seed
set.seed(12345)

#set working directory
setwd("/Users/brycedietrich/Dropbox/workshops/text_as_data/prelim/")

#load nyt_articles
nyt_data<- read.csv("nyt_results.csv",as.is=TRUE,header=FALSE)
names(nyt_data)<-c("date","title","byline","text")

#convert to vector
nyt_articles<-nyt_data$text

#conver to corpus
corpus<-VCorpus(VectorSource(nyt_articles))

#dtm
corpus_dtm <- DocumentTermMatrix(corpus, control = list(stemming = TRUE, stopwords = TRUE, removeNumbers = TRUE, removePunctuation = TRUE, tolower = TRUE))

#term tfidf
#The mean term frequency-inverse document frequency (tf-idf) over documents containing this term is used to select the vocabulary. This measure allows to omit terms which have low frequency as well as those occurring in many documents. We only include terms which have a tf-idf value of at least 0.1 which is a bit more than the median and ensures that the very frequent terms are omitted.
term_tfidf<-tapply(corpus_dtm$v/row_sums(corpus_dtm)[corpus_dtm$i], corpus_dtm$j, mean) * log2(nDocs(corpus_dtm)/col_sums(corpus_dtm > 0))

#eliminate words less than median
corpus_dtm <- corpus_dtm[,term_tfidf >= median(term_tfidf)]
corpus_dtm <- corpus_dtm[row_sums(corpus_dtm) > 0,]

#estimate several topic models
topic_numbers<-seq(2,30)
fitted_models<-NULL
for(topic_number in topic_numbers){
  print(topic_number)
  fitted_models<-c(fitted_models,LDA(corpus_dtm, k = topic_number))	
}

## compute harmonic means
#extract log likelihoods
logLiks<-lapply(fitted_models, function(L){as.numeric(logLik(L))})

#harmonic mean
#http://pages.ucsd.edu/~chl260/PDF/Lee_lda.pdf
harmonicMean <- function(logLikelihoods, precision=2000L) {
  require("Rmpfr")
  llMed <- median(logLikelihoods)
  as.double(llMed - log(mean(exp(-mpfr(logLikelihoods,
                                       prec = precision) + llMed))))
}
h_means<-sapply(logLiks, function(h) harmonicMean(h))

#find number of topics
topic_numbers[which.max(h_means)]

#perplexity
lapply(fitted_models,function(w) perplexity(w))
p_s<-lapply(fitted_models,function(w) perplexity(w))

#find number of topics
topic_numbers[which.max(p_s)]

#entropy
entropy<-sapply(fitted_models, function(x) mean(apply(posterior(x)$topics, 1, function(z) - sum(z * log(z)))))

#find number of topics
topic_numbers[which.max(entropy)]

#let's look at k = 3
#mixtures across documents
posterior(fitted_models[[2]], corpus_dtm)$topics

#distribution of most popular topics
topics(fitted_models[[2]])
table(topics(fitted_models[[2]]))

#five most popular terms
popular_terms <- terms(fitted_models[[2]], 5)