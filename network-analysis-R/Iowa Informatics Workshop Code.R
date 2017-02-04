##################################################################
#Social Networks in R: An Introduction
#Code Prepared for Iowa Informatics Initiative Workshop 
#Instructor: Elizabeth Menninga, University of Iowa
#Thursday, January 12, 2017
##################################################################

##################################################################
#Before the workshop begins, please download the following:
#termDocMatrix.rdata data file from http://www.rdatamining.com/data
#This is Twitter text data from the @rdatamining handle. 
#Thanks to Tanchang Zhao for sharing and making the data available.
###################################################################

###################################################################
#Also, please install the following packages:
#install.packages("igraph") #for analyses and graphics
#install.packages("sna") #for analyses
#install.packages("NetData") #Includes useful datasets for examples
#install.packages("ergm")
#Note: Some of these packages don't work well together, so you'll 
#want to pay attention to which one is loaded. Don't load the packages 
#yet! Just install them.
####################################################################

setwd("C:/Users/emenninga/Dropbox/classes/Networks/I3 Workshop")
load("termDocMatrix.rdata")

# check out part of the matrix
termDocMatrix[5:10,1:20]
#What do we notice about this data right off the bat?

#NOTE: The data is already saved as an R object, but it need not be. You can read in csv as well as excel files easily

# dichotomizing it (don't always need to do this)
termDocMatrix[termDocMatrix>=1] <- 1
termDocMatrix[5:10, 1:20]
#notice the difference in that fifth cell

#Collapsing the two-mode network to a one-mode network. 

#Focusing on terms where relationship is the number of documents the terms have in common 
termMatrix <- termDocMatrix %*% t(termDocMatrix)
# inspect terms numbered 5 to 10
termMatrix[5:10,5:10]

#If we wanted to collapse so that the relationship was number of 
#terms two documents had in common,how could we do that?
    
#########Making network objects in R. 
library(igraph)

# IMPORTANT NOTE: Typically, R objects are numbered
# starting with 1, so if you want the first element in a
# vector, you would reference it by vector_name[1]. HOWEVER,
# igraph objects are numbered starting from 0. 

# build a graph from the above matrix
g <- graph.adjacency(termMatrix, weighted=T, mode = "undirected") 
##Our first real R command. Let's see what it does:
?graph.adjacency

E(g)
V(g)

#If multiple edge attributes can use E(g)$type to look at individual edge attributes.



#Alternate option
g2 <- graph.data.frame(termMatrix, directed=F) 
#Doesn't let you specify weighted, but lets you include vertex "metadata"
#lots of ways to do the same thing in R. It's a blessing and a curse.


######################Vizualizing Network Objects

#Set labels and add an attribute (degree) to vertices
#V tells it to attribute this to the Vertex. 
#E would tell it edges
V(g)$label <- V(g)$name
V(g)$degree <- degree(g)
          
          
# set seed to make the layout reproducible
set.seed(11235) #Why do we need this?

layout1 <- layout.fruchterman.reingold(g)  
##Appropriate for unweighted, undirected graphs.
#Attractive forces occur between adjacent vertices only, 
#whereas repulsive forces occur between every pair of vertices. 
#Each iteration computes the sum of the forces on each vertex, 
#then moves the vertices to their new positions.

plot(g, layout=layout1)

?igraph::layout #For more info on this subcommand. 

g <- simplify(g, remove.loops=T, remove.multiple=F)

plot(g, layout=layout.kamada.kawai)
##kamada.kawai algorithm requires the graph to be connected (only one component and no isolates) 
#Uses distance between nodes to determine where they should go. 
#The shorter the distance between nodes the closer they are.

tkplot(g, layout=layout.kamada.kawai)
#starts with this layout but then you can play around.
          
################## A fancier graph. In R you can control EVERYTHING          
          
V(g)$label.cex <- 2.2 * (V(g)$degree / max(V(g)$degree)) #Scaling the vertex label to reflect the degree
V(g)$label.color <- rgb(0, 0, 1, .8) #very controlling over the color
V(g)$frame.color <- NA #Getting rid of the circle around the nodes
egam <- (log(E(g)$weight)+.4) / max(log(E(g)$weight)+.4) #Scaling the edge weights
E(g)$color <- rgb(1, 0, 0, egam) #color of edges as a function of relative weight
E(g)$width <- egam #width as a function of relative weight

# plot the graph in layout1
plot(g, layout=layout1)



#################Let's do some analysis!
library(NetData) #using this only to get some easy-to-use data.

data(kracknets, package="NetData")

View(krack_full_data_frame)
#Note this is an edge list.

krack_full <- graph.data.frame(krack_full_data_frame) 
summary(krack_full)

#advice relation
E(krack_full)$advice_tie


windows()
plot(krack_full)

tie_type_colors = c(rgb(1,0,0,.5), rgb(0,0,1,.5), rgb(0,0,0,.5))
E(krack_full)$color[ E(krack_full)$advice_tie==1 ] = tie_type_colors[1]
E(krack_full)$color[ E(krack_full)$friendship_tie==1 ] = tie_type_colors[2]
E(krack_full)$color[ E(krack_full)$reports_to_tie==1 ] = tie_type_colors[3]
E(krack_full)$arrow.size=.5 

windows()
plot(krack_full, 
     vertex.color="red", 
     vertex.label=NA, 
     edge.arrow.size=.5)

legend(.5, 1.25,
       legend = c('Advice', 
                  'Friendship',
                  'Reports To'), 
       col = tie_type_colors, 
       lty=1,
       cex = 1)


##########Degree
advice_in <- degree(graph_from_data_frame(advice_data_frame), mode="in") 
advice_out <- degree(graph_from_data_frame(advice_data_frame), mode="out") 
advice_in
advice_out
#What do we notice about this? No, everyone is not connected to everyone.


advice <- subset(advice_data_frame, (advice_tie > 0))
advice_in <- degree(graph.data.frame(advice),mode="in") 
advice_out <- degree(graph.data.frame(advice), mode="out") 
advice_in
advice_out
#Downside to this approach is it removes isolates.


advice_all <- graph.data.frame(advice_data_frame)
advice_thin <- delete.edges(advice_all, E(advice_all)[get.edge.attribute(advice_all,name = "advice_tie")==0])
advice_in1 <- degree(advice_thin, mode="in")
advice_out1 <- degree(advice_thin, mode="out")
advice_in1
advice_out1
##This "removes" edges that are 0 but keeps all the vertices

krack_advice_symmetrized <- as.undirected(advice_thin, mode='collapse')
summary(krack_advice_symmetrized)


friendship <- graph.data.frame(subset(friendship_data_frame, (friendship_tie > 0)))
friendship_in <- degree(friendship, mode="in") 
friendship_out <- degree(friendship, mode="out") 
friendship_in
friendship_out


# we can assemble these measures and export them as a CSV to our working directory.
node_stats_deg <- cbind(advice_in1,
                       advice_out1)

node_stats_deg

write.csv(node_stats_deg, 'krack_node_stats.csv')
##Note: You could export advice_in1 and friendship_in in the same file as well, but would
#need to make sure the nodes are all listed in the same order! They aren't here.
#Sort the data before cbinding and exporting.


######Graph level stuff:
# Density 
graph.density(advice_thin)
graph.density(friendship)

# Reciprocity
reciprocity(advice_thin)
reciprocity(friendship)

# Transitivity
transitivity(advice_thin)
transitivity(friendship)

#Triad Census
triad.census(friendship)


#####Centrality (Node Level)
betweenness_advice <- betweenness(friendship, directed=T) #This is the default, but just to point it out
betweenness_advice

betweenness_advice2 <- betweenness(friendship, directed=F) 
betweenness_advice2 #notice the difference, it symmetrized the matrix

betweenness_advice3 <- betweenness(friendship, directed=T, normalized=T) 
betweenness_advice3


evcent(friendship) #EV Centrality 




#####Community Detection
comm_eb <- edge.betweenness.community(friendship, directed=T) 
comm_eb
comm_eb$bridges

?edge.betweenness.community #Notice it can incorporate weights and assumes directedness
##Note that it automatically maximizes modularity.

plot(as.dendrogram(comm_eb))

###For use later
term_eb<- edge.betweenness.community(g, directed=F, weights=E(g)) #using the Twitter data
plot(as.dendrogram(term_eb))



#######Is this density different from what would be expected if the ties were randomly
#assigned? Are these communities a good partition?

####Need to switch to sna to answer these questions

library(sna)
###NOTE THE ERROR MESSAGES
detach(package:igraph) #Typically I would do this first but I wanted you to see the errors

?cug.test


trans<-cug.test(friendship, gtrans, mode="digraph", cmode="size", reps=1000)
####Note error message. Igraph and sna want the data saved in formats specific to them.

####Let's go back to the twitter data. Object: termMatrix
termMatrix <- termDocMatrix %*% t(termDocMatrix)

###Dichotomizing for ease today
termMatrix[termMatrix>=1] <- 1

dens<-cug.test(termMatrix, gden, mode="graph", cmode="size", reps=1000)
plot.cug.test(dens)

dens_edge<-cug.test(termMatrix, gden, mode="graph", cmode="edges", reps=1000)
dens_edge
plot.cug.test(dens_edge)
##Why does this output make sense?

trans<-cug.test(termMatrix, gtrans, mode="graph", cmode="size", reps=1000)
plot.cug.test(trans)




####################Blockmodeling
# Now we'll look at some community statistics 

# First, we'll set a vector of communities, ordered according  to the order of vertices. 
term_eb$membership 

blockmodel(termMatrix, term_eb$membership)

gden(termMatrix) #for reference

?blockmodel

###NOTE THAT WE HAVE TO SUPPLY THE BLOCKS. This command doesn't identify them for you, but you can give it an hclust object that is
#hierarchical clustering similarity scores and then tell it how many classes to form. Or the height to cut the classes
#But you have to somehow tell it how to determine the blocks.




#################### QAP Examples
#Krackhardt data of connections in high-tech industry

reports<-read.csv("reports.csv")
friends<-read.csv("friends.csv")
advice<-read.csv("advice.csv")
View(reports)


reports<-reports[, -1]
friends<-friends[,-1]
advice<-advice[,-1]


?qaptest
corrtest<-qaptest(list(advice, friends), gcor, g1=1, g2=2, reps=1000)
corrtest
summary(corrtest)


names(corrtest)
corrtest$dist


plot.qaptest(corrtest)


corrtest2<-qaptest(list(reports, friends), gcor, g1=1, g2=2, reps=1000)
corrtest2
summary(corrtest2)


plot.qaptest(corrtest2)




#################Network Regression
?netlogit

classic <- netlogit(advice, list(friends, reports), nullhyp="classical", reps=1000)
summary(classic)

qap <- netlogit(advice, list(friends, reports), nullhyp="qap", reps=1000)
summary(qap)


resultsALL <- cbind(summary(qap)[[1]],summary(classic)[[1]], summary(qap)[[21]], summary(classic)[[20]])
resultsALL






############################Future References
#Don't forget the power of ?command
#https://sna.stanford.edu/rlabs.php includes some very well documented, easy-to-follow code
#The documentation for sna, statnet, and igraph are good places to look to understand all the options within these packages
#And the Journal of Statistical Software includes descriptions of network packages with some examples.
#Good luck!