source('lab02/rand.r')
source('lab02/rtests.r')

rand <- make.rand.lsfr()

sequence <- rand(1000)

monobit.test(sequence)
rev.test(sequence)
