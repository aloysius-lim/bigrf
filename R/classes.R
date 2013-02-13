setClassUnion("character.or.NULL", c("character", "NULL"))
setClassUnion("integer.or.NULL", c("integer", "NULL"))
setClassUnion("numeric.or.NULL", c("numeric", "NULL"))
setClassUnion("logical.or.NULL", c("logical", "NULL"))
setClassUnion("matrix.or.NULL", c("matrix", "NULL"))
setClassUnion("table.or.NULL", c("table", "NULL"))
setClassUnion("big.matrix.or.NULL", c("big.matrix", "NULL"))

setClass("bigcforest",
         representation=representation(
             nexamples="integer",
             varselect="integer",
             factorvars="logical",
             varnlevels="integer",
             contvarseq="integer",
             y="factor",
             ytable="table",
             yclasswts="matrix",
             ntrees="integer",
             nsplitvar="integer",
             maxndsize="integer",
             maxeslevels="integer",
             nrandsplit="integer",
             oobtimes="integer",
             oobvotes="matrix",
             oobpred="integer",
             trainclserr="matrix",
             trainerr="numeric",
             trainconfusion="table",
             varginidec="numeric",
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
             ntrees="integer",
             testytable="table.or.NULL",
             testvotes="matrix",
             testclserr="numeric.or.NULL",
             testerr="numeric.or.NULL",
             testconfusion="table.or.NULL"
         ),
         contains="integer")

setClass("bigrfprox",
         representation=representation(
             examples="big.matrix.or.NULL",
             cachepath="character.or.NULL"
         ),
         contains="big.matrix")
