# Collapsed tree navigation

# The parent of a node in the collapsed tree

#' @keywords internal
#' @export get.tt.parent

get.tt.parent <- function(tt, label){
  if(length(which(tt$unique.split==label))>0){
    tt$parent.split[which(tt$unique.split==label)]
  } else {
    NULL
  }
}

# All ancestors of a node in the collapsed tree.

#' @keywords internal
#' @export get.tt.ancestors

get.tt.ancestors <- function(tt, label){
  out <- vector()
  current.node <- get.tt.parent(tt, label)
  while(current.node != "root"){
    out <- c(out, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  out
}

# All children of a node in the collapsed tree

#' @keywords internal
#' @export get.tt.children

get.tt.children <- function(tt, label){
  tt$unique.split[which(tt$parent.split==label)]
}

# All adjacent nodes of a node in the collapsed tree (parent and children)

#' @keywords internal
#' @export get.tt.adjacent

get.tt.adjacent <- function(tt, label){
  c(get.tt.children(tt, label), get.tt.parent(tt, label))
}

# The grandparent of a node in the collapsed tree

#' @keywords internal
#' @export get.tt.parent.host

get.tt.parent.host <- function(tt, label){
  if(length(which(tt$unique.split==label))>0){
    tt$host[which(tt$unique.split==label)]
  } else {
    NULL
  }
}

# MRCA of a pair of nodes in the collapsed tree

#' @keywords internal
#' @export get.tt.mrca

get.tt.mrca <- function(tt, label1, label2){
  # sanity check
  ca.vec <- intersect(c(label1, get.tt.ancestors(tt, label1)), c(label2, get.tt.ancestors(tt, label2)))
  if(length(ca.vec)>1){
    for(i in seq(1, length(ca.vec)-1)){
      if(get.tt.parent(tt, ca.vec[i]) != ca.vec[i+1]){
        stop("get.tt.mrca not working as intended")
      }
    }
  }
  
  ca.vec[1]
}

# The host corresponding to the MRCA of two nodes

#' @keywords internal
#' @export get.tt.mrca.host

get.tt.mrca.host <- function(tt, label1, label2){
  node <- intersect(c(label1,get.tt.ancestors(tt, label1)), c(label2, get.tt.ancestors(tt, label2)))[1]
  
  tt$host[which(tt$unique.split==node)]
}

# The path from one node in the collapsed tree to another

#' @keywords internal
#' @export get.tt.path

get.tt.path <- function(tt, label1, label2){
  mrca <- get.tt.mrca(tt, label1, label2)
  
  if(mrca == label1){
    result <- vector()
    current.node <- label2
    while(current.node != label1){
      result <- c(result, current.node)
      current.node <- get.tt.parent(tt, current.node)
    }
    result <- c(result, label1)
    return(result)
  }
  
  if(mrca == label2){
    result <- vector()
    current.node <- label1
    while(current.node != label2){
      result <- c(result, current.node)
      current.node <- get.tt.parent(tt, current.node)
    }
    result <- c(result, label2)
    return(result)
  }
  
  first.half <- vector()
  current.node <- label1
  while(current.node != mrca){
    first.half <- c(first.half, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  
  second.half <- vector()
  current.node <- label2
  while(current.node != mrca){
    second.half <- c(second.half, current.node)
    current.node <- get.tt.parent(tt, current.node)
  }
  
  c(first.half, mrca, rev(second.half))
}

# Output the collapsed tree from a phylogeny and its node associations

#' @keywords internal
#' @export output.trans.tree

output.trans.tree <- function(tree, assocs){
  # find the association of each node
  
  assocs.vec <- vector(mode = "character", length = tree$Nnode + length(tree$tip.label))
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(node.no > length(assocs)){
      # annoying and hacky
      assocs.vec[node.no] <- "none"
    } else if (is.null(assocs[[node.no]])){
      assocs.vec[node.no] <- "none"
    } else if (is.na(assocs[[node.no]])){
      assocs.vec[node.no] <- "none"
    } else if (assocs[[node.no]] %in% c("*", "unassigned")) {
      assocs.vec[node.no] <- "none"
    } else {
      if(length(assocs[[node.no]])>1){
        stop("This method is intended for use on trees with conflicts resolved")
      }
      assocs.vec[node.no] <- assocs[[node.no]]
    }
  }
  
  first.of.split <- vector(mode = "logical", length = tree$Nnode + length(tree$tip.label))
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    parent.no <- Ancestors(tree, node.no, type="parent")
    if(parent.no == 0){
      first.of.split[node.no] <- T
    } else {
      first.of.split[node.no] <- assocs.vec[parent.no] != assocs.vec[node.no]
    }
  }
  
  splits.vec <- assocs.vec
  
  unassigned.roots <- which(splits.vec=="none" & first.of.split)
  
  splits.vec[unassigned.roots] <- 
    paste("unassigned_region-SPLIT", 1:length(unassigned.roots), sep="")
  
  for(node.no in seq(1, tree$Nnode + length(tree$tip.label))){
    if(assocs.vec[node.no]=="none" & !first.of.split[node.no]){
      current.node.no <- node.no
      while(!first.of.split[current.node.no]){
        current.node.no <- Ancestors(tree, current.node.no, type="parent")
      }
      if(!grepl("^unassigned_region", splits.vec[current.node.no])){      
        stop("Parent of unassigned node is not one of the unassigned subgraph roots")
      }
      splits.vec[node.no] <- splits.vec[current.node.no]
    }
  }
  
  unique.split <- unique(splits.vec)
  
  tt.table <- tibble(unique.split)
  
  tt.table <- tt.table %>% 
    mutate(root.no = map_int(unique.split, function(x){
      which(splits.vec == x & first.of.split)}), 
      parent.node = Ancestors(tree, root.no, type="parent"),
      parent.split = map_chr(parent.node, function(x){
        if(x == 0){
          "root"
        } else {
          splits.vec[x]
        }
      }),
      length = map2_dbl(root.no, parent.node, function(x, y){
        if(y==0){
          0
        } else {
          tree$edge.length[which(tree$edge[,2]==x)]
        }
      }),
      host = unlist(lapply(strsplit(unique.split, "-SPLIT"), `[[`, 1)),
      parent.host = unlist(lapply(strsplit(parent.split, "-SPLIT"), `[[`, 1))) %>%
    select(unique.split, parent.split, host, parent.host, length, root.no)
  
  tt.table
}

# tt.table %>% mutate(length = map2_chr(root.no, parent.node, function(x, y){
#   if(y==0){
#     0
#   } else {
#     tree$edge.length[which(tree$edge[,2]==x)]
#   }
# }))

#' @keywords internal
#' @export prune.unassigned.tips

prune.unassigned.tips <- function(tt.table){
  
  for.output <- tt.table[,1:6]
  
  unassigned.tips <- which(grepl("unassigned", for.output$unique.splits) &
                             !(for.output$unique.splits %in% for.output$parent.splits))
  
  if(length(unassigned.tips) > 0){
    for.output <- for.output[-unassigned.tips,]
  }
  
  #renumber
  unassigned.rows <- which(grepl("^unassigned_region",for.output$unique.splits))
  
  unassigned.labels <- for.output$unique.splits[which(grepl("^unassigned_region",for.output$unique.splits))]
  
  for(x in 1:length(unassigned.rows)) {
    old.label <- unassigned.labels[x]
    new.label <- paste("UnassignedRegion-SPLIT",x,sep="")
    for.output$unique.splits[unassigned.rows[x]] <- new.label
    for.output$parent.splits[which(for.output$parent.splits==old.label)] <- new.label
  } 
  
  
  for.output$hosts[which(grepl("^unassigned_region",for.output$hosts))] <- "UnassignedRegion"
  
  for.output$parent.hosts[which(grepl("^unassigned_region",for.output$parent.hosts))] <- "UnassignedRegion"
  
  for.output
}

#' @keywords internal
#' @export check.contiguous

check.contiguous <- function(tt, hosts, splits.for.hosts, hosts.for.splits){
  if(length(hosts)!=2){
    stop("Not implemented")
  }
  pat.1.id <- hosts[1]
  pat.2.id <- hosts[2]
  
  OK <- TRUE
  all.nodes <-  c(splits.for.hosts[[pat.1.id]], splits.for.hosts[[pat.2.id]])
  
  for(node.1 in seq(1, length(all.nodes))){
    for(node.2 in seq(1, length(all.nodes))){
      if(node.1 < node.2){
        node.1.id <- all.nodes[node.1]
        node.2.id <- all.nodes[node.2]
        path <- get.tt.path(tt, node.1.id, node.2.id)
        for(node in path){
          if(!grepl("^unassigned_region",node)){
            if(!(hosts.for.splits[[node]] %in% c(pat.1.id, pat.2.id))){
              OK <- FALSE
              break
            }
          }
        }
      }
      if(!OK){
        break
      }
    }
    if(!OK){
      break
    }
  }
  
  OK
}

#' @keywords internal
#' @export check.uninterrupted

check.uninterrupted <- function(tt, hosts, splits.for.hosts, hosts.for.splits){
  # this could certainly be faster
  if(length(hosts)!=2){
    stop("Not implemented")
  }
  pat.1.id <- hosts[1]
  pat.2.id <- hosts[2]
  
  nodes.1 <- splits.for.hosts[[pat.1.id]]
  nodes.2 <- splits.for.hosts[[pat.2.id]]
  
  any.contiguity <- F
  any.interruption <- F
  
  for(node.1 in nodes.1){
    current.node <- get.tt.parent(tt, node.1)
    while(current.node!="root" & (grepl("^unassigned_region",current.node) | (if(is.null(hosts.for.splits[[current.node]])) {T} else {hosts.for.splits[[current.node]]==pat.1.id}))){
      
      current.node <- get.tt.parent(tt, current.node)
    }
    if(current.node != "root"){
      if(hosts.for.splits[[current.node]]==pat.2.id){
        any.contiguity <- T
      } else {
        # this is a blocking node
        chain <- get.tt.ancestors(tt, current.node)
        if(length(intersect(c(pat.1.id, pat.2.id), unlist(hosts.for.splits[chain])))>0){
          any.interruption <- T
          break
        }
      }
    }
  }
  
  if(!any.interruption){
    for(node.2 in nodes.2){
      current.node <- get.tt.parent(tt, node.2)
      while(current.node!="root" & (grepl("unassigned_region",current.node) | (if(is.null(hosts.for.splits[[current.node]])) {T} else {hosts.for.splits[[current.node]]==pat.2.id}))){
        current.node <- get.tt.parent(tt, current.node)
      }
      if(current.node != "root"){
        if(hosts.for.splits[[current.node]]==pat.1.id){
          any.contiguity <- T
        } else {
          # this is a blocking node
          chain <- get.tt.ancestors(tt, current.node)
          if(length(intersect(c(pat.1.id, pat.2.id), unlist(hosts.for.splits[chain])))>0){
            any.interruption <- T
            break
          }
        }
      }
      
    }
  }
  any.contiguity & !any.interruption
}

#' @keywords internal
#' @export extract.tt.subgraph

extract.tt.subgraph <- function(tt, hosts, splits.for.hosts, hosts.for.splits){
  
  if(length(hosts)!=2){
    # for now, at least
    stop("Not implemented")
  }
  if(!check.contiguous(tt, hosts, splits.for.hosts, hosts.for.splits)){
    stop("Not contiguous")
  }
  
  pat.1.id <- hosts[1]
  pat.2.id <- hosts[2]
  
  pat.1.splts <- splits.for.hosts[[pat.1.id]]
  pat.2.splts <- splits.for.hosts[[pat.2.id]]
  
  sub.tt <- tt %>% filter(unique.split %in% c(pat.1.splts, pat.2.splts))
  unassigned.below <- tt[which(tt$host == "unassigned_region" & (tt$parent.split %in% c(pat.1.splts, pat.2.splts))),]
  unassigned.above <- tt[which(tt$host == "unassigned_region" & (tt$unique.split %in% sub.tt$parent.split)),]
  
  none.but.maybe.relevant <- c(unassigned.above$unique.splits, unassigned.below$unique.splits)
  
  adjacent.relevance.count <- sapply(none.but.maybe.relevant, function(x) length(intersect(c(pat.1.splts, pat.2.splts), get.tt.adjacent(tt, x) ) ))
  
  c(sub.tt$unique.split, unique(none.but.maybe.relevant[which(adjacent.relevance.count > 1)]))
}

# every distance between subgraphs

#' @keywords internal
#' @export all.subgraph.distances
#' @importFrom ape dist.nodes node.depth.edgelength

all.subgraph.distances <- function(tree, tt, splits, assocs, slow=F, total.pairs, verbose = F, no.progress.bars = F){
  
  if (verbose & !no.progress.bars) progress.bar <- txtProgressBar(width=50, style=3)
  
  if(!slow){
    tree.dist <- dist.nodes(tree)
  } else {
    depths <- node.depth.edgelength(tree)
  }
  
  count <- 0
  
  out <- matrix(ncol = length(splits), nrow=length(splits))
  
  if(total.pairs==0 & verbose & !no.progress.bars){
    setTxtProgressBar(progress.bar, 1)
  } else {
    for(spt.1.no in 1:length(splits)){
      for(spt.2.no in 1:length(splits)){
        if(spt.1.no==spt.2.no){
          out[spt.1.no, spt.2.no] <- 0
        } else if(spt.1.no<spt.2.no){
          
          
          spt.1 <- splits[spt.1.no]
          spt.2 <- splits[spt.2.no]
          
          chain.1 <- get.tt.ancestors(tt, spt.1)
          chain.2 <- get.tt.ancestors(tt, spt.2)
          
          if(spt.1 %in% chain.2){
            mrca.2 <- tt$root.no[which(tt$unique.split==spt.2)]
            current.node <- mrca.2
            length <- 0
            while(assocs[[current.node]]!=spt.1){
              length <- length + get.edge.length(tree, current.node)
              current.node <- Ancestors(tree, current.node, type="parent")
              if(is.root(tree, current.node)){
                stop("Reached the root?")
              }
            }
            out[spt.1.no, spt.2.no] <- length
            
          } else if(spt.2 %in% chain.1){
            mrca.1 <- tt$root.no[which(tt$unique.split==spt.1)]
            current.node <- mrca.1
            length <- 0
            while(assocs[[current.node]]!=spt.2){
              
              length <- length + get.edge.length(tree, current.node)
              current.node <- Ancestors(tree, current.node, type="parent")
              if(is.root(tree, current.node)){
                stop("Reached the root?")
              }
            }
            out[spt.1.no, spt.2.no] <- length
            
          } else {
            
            mrca.1 <- tt$root.no[which(tt$unique.split==spt.1)]
            mrca.2 <- tt$root.no[which(tt$unique.split==spt.2)]
            if(!slow){
              out[spt.1.no, spt.2.no] <- tree.dist[mrca.1, mrca.2]
            } else {
              out[spt.1.no, spt.2.no] <- pat.dist(tree, depths, mrca.1, mrca.2)
            }
          }
          
          count <- count + 1
          
          if(verbose & !no.progress.bars){
            setTxtProgressBar(progress.bar, count/total.pairs)
          }
          out[spt.2.no, spt.1.no] <- out[spt.1.no, spt.2.no]
        }
      }
    }
  }
  
  if(verbose & !no.progress.bars) close(progress.bar)
  colnames(out) <- splits
  rownames(out) <- splits
  
  out
}

#' @keywords internal
#' @export pat.dist

pat.dist <- function(tree, depths, node.1, node.2){
  node.1.chain <- Ancestors(tree, node.1, type="all")
  node.2.chain <- Ancestors(tree, node.2, type="all")
  if(node.1 %in% node.2.chain){
    mrca <- node.1
    dist <- depths[node.2] - depths[node.1]
  } else if(node.2 %in% node.1.chain){
    mrca <- node.2
    dist <- depths[node.1] - depths[node.2]
  } else {
    mrca <- node.1.chain[which(node.1.chain %in% node.2.chain)][1]
    dist <- (depths[node.1] - depths[mrca]) + (depths[node.2] - depths[mrca])
  }
  
  dist
}

# are each pair of subgraphs adjacent? If !none.matters then two nodes separated only by "none" are still adjacent

#' @keywords internal
#' @export subgraphs.adjacent

subgraphs.adjacent <- function(tt, splits, none.matters = F){
  out <- matrix(ncol = length(splits), nrow=length(splits))
  for(spt.1.no in 1:length(splits)){
    for(spt.2.no in 1:length(splits)){
      if(spt.1.no==spt.2.no){
        out[spt.1.no, spt.2.no] <- NA
      } else if(spt.1.no<spt.2.no){
        spt.1 <- splits[spt.1.no]
        spt.2 <- splits[spt.2.no]
        
        if(spt.1 %in%  get.tt.adjacent(tt, spt.2)){
          out[spt.1.no, spt.2.no] <- T
          out[spt.2.no, spt.1.no] <- T
        } else if(!none.matters) {
          path <- get.tt.path(tt, spt.1, spt.2)
          internal.path <- path[2:(length(path)-1)]
          adj <- length(internal.path)==1 & grepl("unassigned_region",internal.path[1])
          out[spt.1.no, spt.2.no] <- adj
          out[spt.2.no, spt.1.no] <- adj
        } else {
          out[spt.1.no, spt.2.no] <- F
          out[spt.2.no, spt.1.no] <- F
        }
      }
    }
  }
  colnames(out) <- splits
  rownames(out) <- splits
  
  out
}

# are pairs of subgraphs from two hosts not separated by any other subgraphs from either of those hosts?

#' @keywords internal
#' @export subgraphs.unblocked

subgraphs.unblocked <- function(tt, splits, total.pairs, verbose = F, no.progress.bars = F){
  
  if (verbose & !no.progress.bars) progress.bar <- txtProgressBar(width=50, style=3)
  count <- 0
  
  out <- matrix(ncol = length(splits), nrow=length(splits))
  if(total.pairs==0 & verbose & !no.progress.bars){
    setTxtProgressBar(progress.bar, 1)
  } else {
    
    for(spt.1.no in 1:length(splits)){
      for(spt.2.no in 1:length(splits)){
        if(spt.1.no==spt.2.no){
          out[spt.1.no, spt.2.no] <- NA
        } else if(spt.1.no<spt.2.no){
          spt.1 <- splits[spt.1.no]
          spt.2 <- splits[spt.2.no]
          
          pat.1 <- strsplit(spt.1, "-SPLIT")[[1]][1]
          pat.2 <- strsplit(spt.2, "-SPLIT")[[1]][1]
          
          if(spt.1 %in% get.tt.adjacent(tt, spt.2)){
            out[spt.1.no, spt.2.no] <- T
            out[spt.2.no, spt.1.no] <- T
          } else {
            path <- get.tt.path(tt, spt.1, spt.2)
            internal.path <- path[2:(length(path)-1)]
            blockers <- which(grepl(paste0('^',pat.1),internal.path) | grepl(paste0('^',pat.2),internal.path))
            
            adj <- length(blockers) == 0
            out[spt.1.no, spt.2.no] <- adj
            out[spt.2.no, spt.1.no] <- adj
          } 
          count <- count + 1
          if (verbose & !no.progress.bars) {
            setTxtProgressBar(progress.bar, count/total.pairs)
          }
        }
      }
    }
  }
  if (verbose & !no.progress.bars) close(progress.bar)
  
  colnames(out) <- splits
  rownames(out) <- splits
  
  out
}


#' @keywords internal
#' @export check.adjacency

check.adjacency <- function(tt, hosts, splits.for.hosts){
  # this could certainly be faster
  if(length(hosts)!=2){
    stop("Not implemented")
  }
  pat.1.id <- hosts[1]
  pat.2.id <- hosts[2]
  
  nodes.1 <- splits.for.hosts[[pat.1.id]]
  nodes.2 <- splits.for.hosts[[pat.2.id]]
  
  for(node.1 in nodes.1){
    for(node.2 in nodes.2){
      if(check.tt.node.adjacency(tt, node.1, node.2, T)){
        return(T)
      }
    }
  } 
  
  F
}

#' @keywords internal
#' @export check.tt.node.adjacency

check.tt.node.adjacency <- function(tt, label1, label2, allow.unassigned = F){
  path <- get.tt.path(tt, label1, label2)
  
  if(!allow.unassigned){
    return(length(path)==2)
  }
  
  if(length(path)==2){
    # they must be next to each other
    return(T)
  }
  
  if(length(path)>3){
    # they can't be - a path length greater than 2 can only go through an unassigned region, and adjacent unassigned regions are not allowed
    return(F)
  }
  
  grepl("^unassigned_region",path[2])
}

#' @keywords internal
#' @export classify

classify <- function(ptree, allow.mt = F, n.mt=Inf, p.mt= Inf, zero.length.adjustment=F, relaxed.ancestry = F,verbose = F, no.progress.bars = F) {	
  
  if(is.null(ptree[["tree"]])){
    
    if (verbose) cat("Reading tree file ", ptree$tree.file.name, "...\n", sep = "")
    
    pseudo.beast.import <- read.beast(ptree$tree.file.name)
    
    if (verbose) cat("Reading annotations...\n")
    
    tree <- as.phylo(pseudo.beast.import)
    
    annotations <- pseudo.beast.import[,c("node", "INDIVIDUAL", "SPLIT")]
  
    # to deal with a problem with older versions of ggtree
    
    strip.quotes <- function(string) if(substr(string, 1, 1)=="\"") substr(string, 2, nchar(string)-1) else string
    
    annotations$INDIVIDUAL <- sapply(annotations$INDIVIDUAL, strip.quotes)
    annotations$SPLIT <- sapply(annotations$SPLIT, strip.quotes)
    
  } else {
    
    tree <- ptree$tree
    
    annotations <- tibble(node = seq(1, length(tree$tip.label) + tree$Nnode), 
                              INDIVIDUAL = as.character(attr(tree, "INDIVIDUAL")), 
                              SPLIT = as.character(attr(tree, "SPLIT")))
  }
  
  
  if(is.null(ptree$splits.table)){
    if (verbose) cat("Reading splits file", ptree$splits.file.name, "...\n")
    
    splits <- read_csv(ptree$splits.file.name)
  } else {
    splits <- ptree$splits.table
  }
  
  if (verbose) cat("Collecting tips for each host...\n")
  
  hosts <- unique(splits$host)
  
  hosts <- hosts[hosts!="unassigned"]
  
  all.splits <- unique(splits$subgraph)
  all.splits <- all.splits[all.splits!="unassigned"]
  
  in.order <- match(seq(1, length(tree$tip.label) + tree$Nnode), annotations$node)
  
  assocs <- annotations$SPLIT[in.order]
  
  assocs <- lapply(assocs, function(x) replace(x, is.na(x), "none"))
  assocs <- lapply(assocs, function(x) replace(x, x=="unassigned", "none"))
  
  splits.for.hosts <- lapply(hosts, function(x) unique(splits$subgraph[which(splits$host==x)] ))
  names(splits.for.hosts) <- hosts
  
  hosts.for.splits <- lapply(all.splits, function(x) unique(splits$host[which(splits$subgraph==x)] ))
  names(hosts.for.splits) <- all.splits
  
  hosts.included <- hosts
  
  total.split.pairs <- (length(all.splits)^2 - length(all.splits))/2
  total.host.pairs <- (length(hosts)^2 - length(hosts))/2
  
  individual <- attr(tree, 'INDIVIDUAL')
  
  if (verbose) cat("Collapsing subgraphs...\n")
  
  tt <- output.trans.tree(tree, assocs)
  
  if (verbose) cat("Identifying pairs of unblocked splits...\n")
  
  collapsed.adjacent <- subgraphs.unblocked(tt, all.splits, total.split.pairs, verbose, no.progress.bars)
  
  if (verbose) cat("Calculating pairwise distances between splits...\n")
  
  split.distances <- tryCatch(
    all.subgraph.distances(tree, tt, all.splits, assocs, FALSE, total.split.pairs, verbose, no.progress.bars), warning=function(w){return(NULL)}, error=function(e){return(NULL)})
  
  if(is.null(split.distances)){
    split.distances <- all.subgraph.distances(tree, tt, all.splits, assocs, TRUE, total.split.pairs, verbose, no.progress.bars)
  }
  
  if (verbose) cat("Testing pairs...\n")
  
  progress.bar <- NULL
  
  if (verbose & !no.progress.bars) progress.bar <- txtProgressBar(width=50, style=3)
  
  count <- 0
  adjacency.matrix <- matrix(NA, length(hosts.included), length(hosts.included), dimnames=list(hosts.included, hosts.included))
  contiguity.matrix <- matrix(NA, length(hosts.included), length(hosts.included), dimnames=list(hosts.included, hosts.included))
  top.class.matrix <- matrix(NA, length(hosts.included), length(hosts.included), dimnames=list(hosts.included, hosts.included))
  nodes.1.matrix <- matrix(NA, length(hosts.included), length(hosts.included), dimnames=list(hosts.included, hosts.included))
  nodes.2.matrix <- matrix(NA, length(hosts.included), length(hosts.included), dimnames=list(hosts.included, hosts.included))
  dir.12.matrix <- matrix(NA, length(hosts.included), length(hosts.included), dimnames=list(hosts.included, hosts.included))
  dir.21.matrix <- matrix(NA, length(hosts.included), length(hosts.included), dimnames=list(hosts.included, hosts.included))
  min.distance.matrix <- matrix(NA, length(hosts.included), length(hosts.included), dimnames=list(hosts.included, hosts.included))
  
  
  n.mt <- 3
  
  if(total.host.pairs==0 & verbose & !no.progress.bars){
    setTxtProgressBar(progress.bar, 1)
  } else {
    
    for(pat.1 in seq(1, length(hosts.included))){
      for(pat.2 in  seq(1, length(hosts.included))){
        if (pat.1 < pat.2) {
          
          count <- count + 1
          
          pat.1.id <- hosts.included[pat.1]
          pat.2.id <- hosts.included[pat.2]
          
          nodes.1 <- splits.for.hosts[[pat.1.id]]
          nodes.2 <- splits.for.hosts[[pat.2.id]]
          
          all.nodes <-  c(nodes.1, nodes.2)
          
          adjacency.matrix[pat.1, pat.2] <- check.adjacency(tt, c(pat.1.id, pat.2.id), splits.for.hosts)		
          contiguity.matrix[pat.1, pat.2] <- check.contiguous(tt, c(pat.1.id, pat.2.id), splits.for.hosts, hosts.for.splits)		
          
          count.12 <- 0
          count.21 <- 0
          
          for(node.2 in nodes.2){
            ancestors <- get.tt.ancestors(tt, node.2)
            if(length(intersect(ancestors, nodes.1)) > 0){
              count.12 <- count.12 + 1
            }
          }
          
          for(node.1 in nodes.1){
            ancestors <- get.tt.ancestors(tt, node.1)
            if(length(intersect(ancestors, nodes.2)) > 0){
              count.21 <- count.21 + 1
            }
          }
          
          prop.12 <- count.12/length(nodes.2)
          prop.21 <- count.21/length(nodes.1)
          
          nodes.1.matrix[pat.1, pat.2] <- length(nodes.1)
          nodes.2.matrix[pat.1, pat.2] <- length(nodes.2)
          
          dir.12.matrix[pat.1, pat.2] <- count.12
          dir.21.matrix[pat.1, pat.2] <- count.21
        
          multifurcation.threshold <- 1e-5
          depths <- node.depth.edgelength(tree)
          

          if(1)
          {   
            # tip counts
            tmp <- grep('_read_[0-9]+_count_[0-9]+$',tree$tip.label,value = T)
            count.tip.1 <- length(grep(pat.1.id,tmp,value=T))
            count.tip.2 <- length(grep(pat.2.id,tmp,value=T))
            p.12 <- count.12/count.tip.2
            p.21 <- count.21/count.tip.1
            
            
            if(n.mt<1e5 & p.mt<=1){
              warning('n.mt and p.mt were all specified, and only p.mt was used.')
            }
            
            if(count.12 == 0 & count.21 == 0){
              top.class.matrix[pat.1, pat.2] <- "noAncestry"
            } else if(count.12 != 0 & count.21 == 0 & (relaxed.ancestry | prop.12 == 1)) {
              if(count.12 == 1){
                top.class.matrix[pat.1, pat.2] <- "anc"
              } else {
                if(allow.mt){
                  # p.mt
                  if(p.mt<=1){
                    if(p.12 >=p.mt){
                      top.class.matrix[pat.1, pat.2] <- "complex"
                    }else{
                      top.class.matrix[pat.1, pat.2] <- "multiAnc"
                    }
                  # n.mt
                  }else if(n.mt<1e5){
                    if(count.12>=n.mt){
                      top.class.matrix[pat.1, pat.2] <- "complex"
                    }else{
                      top.class.matrix[pat.1, pat.2] <- "multiAnc"
                    }
                  # zero length adjustment
                  }else if(zero.length.adjustment==T){
                    anc.nodes <- which(attr(ptree$tree,"INDIVIDUAL")==pat.1.id & attr(ptree$tree, "SUBGRAPH_MRCA")==T)
                    desc.nodes <- which(attr(ptree$tree,"INDIVIDUAL")==pat.2.id & attr(ptree$tree, "SUBGRAPH_MRCA")==T)
                    anc.tips <- ptree$tips.for.hosts[[pat.1.id]]
                    anc.mrca.node <- anc.nodes[!(anc.nodes %in% anc.tips)]
                    cb <- expand.grid(anc.mrca.node,desc.nodes)
                    desc.anc.pd <- sapply(1:nrow(cb),
                                          function(x)pat.dist(tree, depths, cb$Var1[x], cb$Var2[x]))
                    if(any(desc.anc.pd<multifurcation.threshold, na.rm = T)){
                      top.class.matrix[pat.1, pat.2] <- "complex"
                    }else{
                      top.class.matrix[pat.1, pat.2] <- "multiAnc"
                    }
                  }else{
                    top.class.matrix[pat.1, pat.2] <- "multiAnc"
                  }
                } else {
                  top.class.matrix[pat.1, pat.2] <- "complex"
                }
              }
            } else if(count.21 != 0 & count.12 == 0 & (relaxed.ancestry | prop.21 == 1)) {
              if(count.21 == 1){
                top.class.matrix[pat.1, pat.2] <- "desc"
              } else {
                if(allow.mt){
                  # p.mt
                  if(p.mt<=1){
                    if(p.21 >=p.mt){
                      top.class.matrix[pat.1, pat.2] <- "complex"
                    }else{
                      top.class.matrix[pat.1, pat.2] <- "multiDesc"
                    }
                  # n.mt
                  }else if(n.mt<1e5){
                    if(count.21>=n.mt){
                      top.class.matrix[pat.1, pat.2] <- "complex"
                    }else{
                      top.class.matrix[pat.1, pat.2] <- "multiDesc"
                    }
                  # zero length adjustment
                  }else if(zero.length.adjustment==T){
                    anc.nodes <- which(attr(ptree$tree,"INDIVIDUAL")==pat.2.id & attr(ptree$tree, "SUBGRAPH_MRCA")==T)
                    desc.nodes <- which(attr(ptree$tree,"INDIVIDUAL")==pat.1.id & attr(ptree$tree, "SUBGRAPH_MRCA")==T)
                    anc.tips <- ptree$tips.for.hosts[[pat.2.id]]
                    anc.mrca.node <- anc.nodes[!(anc.nodes %in% anc.tips)]
                    cb <- expand.grid(anc.mrca.node,desc.nodes)
                    desc.anc.pd <- sapply(1:nrow(cb),
                                          function(x)pat.dist(tree, depths, cb$Var1[x], cb$Var2[x]))
                    if(any(desc.anc.pd<multifurcation.threshold, na.rm = T)){
                      top.class.matrix[pat.1, pat.2] <- "complex"
                    }else{
                      top.class.matrix[pat.1, pat.2] <- "multiDesc"
                    }
                  }else{
                    top.class.matrix[pat.1, pat.2] <- "multiDesc"
                  }
                } else {
                  top.class.matrix[pat.1, pat.2] <- "complex"
                }
              }
            } else {
              top.class.matrix[pat.1, pat.2] <- "complex"
            }
          }
          
          pairwise.distances <- vector()
          
          for(node.1 in nodes.1){
            for(node.2 in nodes.2){
              if(collapsed.adjacent[node.1, node.2]){
                pairwise.distances <- c(pairwise.distances, split.distances[node.1, node.2])
              }
            }
          }
          
          min.distance.matrix[pat.1, pat.2] <- min(pairwise.distances)
          
          if (verbose & !no.progress.bars) {
            setTxtProgressBar(progress.bar, count/total.host.pairs)
          }
        }
      }
    }
  }
  
  if (verbose & !no.progress.bars) close(progress.bar)
  
  normalisation.constant <- ptree$normalisation.constant
  
  # get the first two columns once only
  # this procedure still turns strings into factors for reasons I don't really understand
  
  melted.matrix <- adjacency.matrix %>% melt() %>% filter(!is.na(value)) 
  
  host.1 <- melted.matrix %>% magrittr::extract(, "Var1") %>% as.character()
  host.2 <- melted.matrix %>% magrittr::extract(, "Var2") %>% as.character()
  adjacent <- melted.matrix %>% magrittr::extract(, "value")
  
  contiguous <- contiguity.matrix %>% melt() %>% filter(!is.na(value)) %>% magrittr::extract(, "value")
  paths12 <- dir.12.matrix %>% melt() %>% filter(!is.na(value)) %>% magrittr::extract(, "value")
  paths21 <- dir.21.matrix %>% melt() %>% filter(!is.na(value)) %>% magrittr::extract(, "value")
  nodes1 <- nodes.1.matrix %>% melt() %>% filter(!is.na(value)) %>% magrittr::extract(, "value")
  nodes2 <- nodes.2.matrix %>% melt() %>% filter(!is.na(value)) %>% magrittr::extract(, "value")
  ancestry <- top.class.matrix %>% melt() %>% filter(!is.na(value)) %>% magrittr::extract(, "value") %>% as.character()
  min.distance.between.subgraphs <- min.distance.matrix %>% melt() %>% filter(!is.na(value)) %>% magrittr::extract(, "value")

  classification <- tibble(host.1, host.2, adjacent, contiguous, paths12, paths21, nodes1, nodes2, ancestry, min.distance.between.subgraphs)
  
  if(normalisation.constant!=1){
    normalised.distance.matrix <- min.distance.matrix/normalisation.constant
    normalised.min.distance.between.subgraphs <- normalised.distance.matrix %>% melt() %>% filter(!is.na(value)) %>% magrittr::extract(, "value")
    classification <- classification %>% mutate(normalised.min.distance.between.subgraphs = normalised.min.distance.between.subgraphs)
  }
  
  return(list(classification = classification, collapsed=tt))
}

#' @keywords internal
#' @export merge.classifications


merge.classifications <- function(ptrees, verbose = F){
  classification.rows	<- ptrees %>% map(function(ptree) {
    
    if(is.null(ptree$classification.results$classification) & is.null(ptree$classification.file.name)){
      NULL
    }
    
    if(!is.null(ptree$classification.results$classification)){
      # TODO this doesn't need to be coerced to a tibble in the end - it should be already
      
      tt <- as_tibble(ptree$classification.results$classification)
      
      # TODO remove this hack
      if("path.classification" %in% names(tt)) tt <- tt %>% rename(ancestry = path.classification)
    } else {
      if (verbose) cat("Reading window input file ", ptree$classification.file.name,"\n", sep="")
      tt <- read_csv(ptree$classification.file.name, col_names = TRUE)
      
      if(nrow(tt)==0){
        stop(paste0("No lines in window input file ",ptree$classification.file.name))
      }
    }
    
    tt <- tt %>% add_column(tree.id = ptree$id)
    
    tt
  })	
  if(length(classification.rows)==0){
    stop("No classification results present in any window; cannot continue.\n")
  }
  
  # need to worry about transmissions listed in the wrong direction in the input file
  
  if(verbose) cat("Rearranging host pairs...\n")
  
  classification.rows	<- classification.rows %>% map(function(x){    
    x %>%
        mutate(tempancestry = ancestry, temphost.1 = host.1, temphost.2 = host.2, temppaths12 = paths12, temppaths21 = paths21, tempnodes1 = nodes1, tempnodes2 = nodes2) %>%
        mutate(ancestry = replace(ancestry, host.2 < host.1 & tempancestry == "anc", "desc")) %>% 
        mutate(ancestry = replace(ancestry, host.2 < host.1 & tempancestry == "desc", "anc")) %>% 
        mutate(ancestry = replace(ancestry, host.2 < host.1 & tempancestry == "multiAnc", "multiDesc")) %>% 
        mutate(ancestry = replace(ancestry, host.2 < host.1 & tempancestry == "multiDesc", "multiAnc")) %>%
        mutate(host.1 = replace(host.1, temphost.2 < temphost.1, temphost.2[temphost.2 < temphost.1])) %>%
        mutate(host.2 = replace(host.2, temphost.2 < temphost.1, temphost.1[temphost.2 < temphost.1])) %>%
        mutate(paths12 = replace(paths12, temphost.2 < temphost.1, temppaths21[temphost.2 < temphost.1])) %>%
        mutate(paths21 = replace(paths21, temphost.2 < temphost.1, temppaths12[temphost.2 < temphost.1])) %>%
        mutate(nodes1 = replace(nodes1, temphost.2 < temphost.1, tempnodes2[temphost.2 < temphost.1])) %>%
        mutate(nodes2 = replace(nodes2, temphost.2 < temphost.1, tempnodes1[temphost.2 < temphost.1])) %>%
        select(-c(temphost.1, temphost.2, tempancestry, temppaths12, temppaths21, tempnodes1, tempnodes2))
  })
  #
  # rbind consolidated files
  #

  if(verbose) cat("Consolidating file contents...\n")

  classification.rows <- bind_rows(classification.rows)

  if(verbose) cat("Finding patristic distance columns...\n")
  
  # reset names depending on which classify script was used
  if('normalised.min.distance.between.subgraphs' %in% names(classification.rows)){
    
    classification.rows <- classification.rows %>% rename(patristic.distance = normalised.min.distance.between.subgraphs)
    classification.rows <- classification.rows %>% select(-min.distance.between.subgraphs)
    
  } else if('min.distance.between.subgraphs' %in% names(classification.rows)){
    
    classification.rows <- classification.rows %>% rename(patristic.distance = min.distance.between.subgraphs)
    
  }
  
  #	check we have patristic distances, paths
  
  stopifnot('patristic.distance' %in% names(classification.rows))
  stopifnot(!nrow(classification.rows %>% filter(is.na(patristic.distance))))
  stopifnot(!nrow(classification.rows %>% filter(is.na(paths12))))
  stopifnot(!nrow(classification.rows %>% filter(is.na(paths21))))
  
  #	set to numeric (probably not necessary now)
  
  classification.rows <- classification.rows %>% mutate_at("patristic.distance", as.numeric)
  
  if(verbose) cat("Reordering...\n")
  
  #	reorder
  classification.rows <- classification.rows %>% arrange(tree.id, host.1, host.2)
  
  return(classification.rows)
}


