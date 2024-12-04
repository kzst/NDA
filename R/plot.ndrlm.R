#-----------------------------------------------------------------------------#
#                                                                             #
#  GENERALIZED NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (GNDA)     #
#                                                                             #
#  Written by: Zsolt T. Kosztyan*, Marcell T. Kurbucz, Attila I. Katona,      #
#              Zahid Khan                                                     #
#              *Department of Quantitative Methods                            #
#              University of Pannonia, Hungary                                #
#              kosztyan.zsolt@gtk.uni-pannon.hu                               #
#                                                                             #
# Last modified: February 2024                                                #
#-----------------------------------------------------------------------------#
###### PLOT FOR NETWORK-BASED DIMENSIONALITY REDUCTION AND ANALYSIS (NDA) #####
#' @export
plot.ndrlm <- function(x,sig=0.05,interactive=FALSE,...){
  if (methods::is(x,"ndrlm")){
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop(
        "Package \"igraph\" must be installed to use this function.",
        call. = FALSE
      )
    }
    if (!requireNamespace("stats", quietly = TRUE)) {
      stop(
        "Package \"stats\" must be installed to use this function.",
        call. = FALSE
      )
    }
    if (!requireNamespace("visNetwork", quietly = TRUE)) {
      stop(
        "Package \"visNetwork\" must be installed to use this function.",
        call. = FALSE
      )
    }
    if (!requireNamespace("lm.beta", quietly = TRUE)) {
      stop(
        "Package \"lm.beta\" must be installed to use this function.",
        call. = FALSE
      )
    }
    nY<-ncol(x$Y)
    nS<-ncol(x$NDA$scores)
    nX<-ncol(x$X)
    membership<-x$NDA$membership
    loadings<-x$NDA$loadings
    node_ID<-1:(nY+nS+nX)
    node_label<-c(colnames(x$Y),colnames(x$NDA$scores),colnames(x$X))
    node_shape<-c(rep("rectangle",nY),rep("circle",nS),
                  rep("rectangle",nX))
    nodes=data.frame(id=node_ID,label=node_label,shape=node_shape)
    nodes_color<-c(rep(0,nY),1:nS,membership)
    nodes<-data.frame(id=node_ID,label=node_label,shape=node_shape,
                      color=nodes_color)
    edges <- data.frame(matrix(ncol = 3, nrow = 0))
    colnames(edges) <- c('from', 'to', 'weight')
    k<-1
    for (i in 1:nY){
      coefs<-as.vector(lm.beta::lm.beta(x$fits[[i]])$standardized.coefficients)[-1]
      pvalues<-summary(x$fits[[i]])$coefficients[-1,4]
      for (j in 1:length(coefs)){
        if (pvalues[j]<sig){
          edges[k,"to"]<-i
          edges[k,"from"]<-nY+j
          edges[k,"weight"]<-coefs[j]
          k<-k+1
        }
      }
    }

    for (i in 1:nX){
      for (j in 1:nS){
        if (membership[i]==j){
          edges[k,"to"]<-nY+j
          edges[k,"from"]<-nY+nS+i
          edges[k,"weight"]<-loadings[colnames(x$X)[i],j]
          k<-k+1
        }
      }
    }
    space=150
    cust_layout<-matrix(0,ncol=2,nrow=nY+nS+nX)
    cust_layout[1:nY,1]<-2
    cust_layout[1:nY,2]<-((1:nY)-mean(1:nY))*space
    cust_layout[(nY+1):(nY+nS),1]<-1
    cust_layout[(nY+1):(nY+nS),2]<-((1:nS)-mean(1:nS))*space
    cust_layout[(nY+nS+1):(nY+nS+nX),1]<-0
    cust_layout[(nY+nS+1):(nY+nS+nX),2]<-((1:nX)-mean(1:nX))*space

    G<-igraph::graph_from_data_frame(edges,
                             directed=TRUE,
                             vertices=nodes)

    if (interactive==TRUE){
      nodes$color<-grDevices::hsv((c(rep(0,nY),1:nS,membership)+1)/
                                    max((c(rep(0,nY),1:nS,membership)+1)),
                                  alpha=0.4)
      edges$arrows=ifelse(igraph::is.directed(G),c("to"),"")
      edges$width=(abs(igraph::E(G)$weight))
      nodes$shape<-gsub("rectangle","box",nodes$shape)
      nodes$shape<-gsub("circle","ellipse",nodes$shape)
      edges$label<-as.vector(paste(round(edges$weight,2),sep=""))
      nw <-
        visNetwork::visIgraphLayout(
          visNetwork::visNodes(
            visNetwork::visInteraction(
              visNetwork::visOptions(
                visNetwork::visEdges(
                  visNetwork::visNetwork(
                    nodes, edges, height = "1000px", width = "100%"),
                  font = list(size = 6),color="#555555",
                  label=edges$label),
                highlightNearest = TRUE, selectedBy = "label"),
              dragNodes = TRUE,
              dragView = TRUE,
              zoomView = TRUE,
              hideEdgesOnDrag = FALSE),physics=FALSE, size=16,
            borderWidth = 1,
            shape=nodes$shape,
            font=list(face="calibri")),layout="layout.norm",
          layoutMatrix = cust_layout,
          physics = FALSE, type="full"
        )
      nw

    }else{
      plot(G,layout=cust_layout,edge.width=abs(igraph::E(G)$weight)*10,
           edge.label=round(igraph::E(G)$weight,2))
    }

  }else{
    plot(x,...)
  }
}
