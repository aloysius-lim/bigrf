setClassUnion("character.or.NULL", c("character", "NULL"))
setClassUnion("integer.or.NULL", c("integer", "NULL"))
setClassUnion("numeric.or.NULL", c("numeric", "NULL"))
setClassUnion("logical.or.NULL", c("logical", "NULL"))
setClassUnion("table.or.NULL", c("table", "NULL"))

setClass("bigcforest",
         representation=representation(
             supervised="logical",
             nexamples="integer",
             factors="logical",
             ylevels="character",
             nlevels="integer",
             ytable="table",
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
             cachepath="character.or.NULL"
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

setClass("bigcprediction",
         representation=representation(
             ntest="integer",
             testlabelled="logical",
             nclass="integer",
             ntrees="integer",
             testytable="table.or.NULL",
             testvotes="matrix",
             testclserr="numeric.or.NULL",
             testerr="numeric.or.NULL",
             testconfusion="table.or.NULL"
         ),
         contains="integer")
