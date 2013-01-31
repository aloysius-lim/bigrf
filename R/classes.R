setClass("bigcforest",
         representation=representation(
           # x
           # y
           # asave
           supervised="logical",
           nsample="integer",
           factors="logical",
           ylevels="character",
           nlevels="integer",
           varselect="integer",
           contvarseq="integer",
           nclass="integer",
           classweights="matrix",
           ntrees="integer",
           maxnodes="integer",
           nsplitvar="integer",
           maxndsize="integer",
           maxeslevels="integer",
           nrandsplit="integer",
           oobtimes="integer",
           oobvotes="matrix",
           oobpred="integer",
           trainclserr="numeric",
           trainerr="numeric",
           trainconfusion="table",
           avgini="numeric",
           printerrfreq="integer",
           printclserr="logical",
           cachepath="character"
         ),
         contains="list")

setClass("bigctree",
         representation=representation(
           insamp="integer",
           inweight="numeric",
           nnodes="integer",
           treemap="matrix",
           nodeclass="integer",
           nodewt="numeric",
           bestvar="integer",
           bestnumsplit="numeric",
           bestcatsplit="list",
           termincount="numeric",
           trainprednode="integer",
           trainpredclass="integer",
           tgini="numeric"
         ))

