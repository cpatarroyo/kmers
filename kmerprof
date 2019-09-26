#Required sources
using FastaIO
using Statistics
#using StatsBase

#Read multi-fasta file
sequences = readfasta("Viruses.fasta")

#Constant definition
counter=1
kmer=4
letters=["A","C","G","T"]
zval=zeros(Float32,(4^kmer,length(sequences)))

#Generate all permutations of length n
if kmer==4
    kmin2=map(x -> join(x),Iterators.product(letters,letters))
elseif kmer==6
    kmin2=map(x -> join(x),Iterators.product(letters,letters,letters,letters))
elseif kmer==8
    kmin2=map(x -> join(x),Iterators.product(letters,letters,letters,letters,letters,letters))
end
kmin1=map(x -> join(x), Iterators.product(kmin2,letters))
kmers=map(x -> join(x), Iterators.product(kmin1,letters))

for sequ in sequences

    #Initialize the dictionaries
    dictn=Dict{String,Int16}()
    dictmin1=Dict{String,Int16}()
    dictmin2=Dict{String,Int16}()
    dictres=Dict{String,Float16}()

    #Count the frequency of the k-mers, k-1mers and k-2mers to calculate z-values
    for k in collect(kmin2)
        dictmin2[k]=length(collect(eachmatch(Regex(k), sequ[2], overlap=true)))
    end
    for j in collect(kmin1)
        dictmin1[j]=length(collect(eachmatch(Regex(j), sequ[2], overlap=true)))
    end
    for i in collect(kmers)
        dictn[i]=length(collect(eachmatch(Regex(i), sequ[2], overlap=true)))
        #Calculate the z-value for each k-mer
        if dictmin2[SubString(i,2,1+(kmer-2))] == 0
            dictres[i] = 0
        else
            dictres[i] = (dictn[i] - ((dictmin1[SubString(i,1,(kmer-1))]*dictmin1[SubString(i,2,kmer)])/dictmin2[SubString(i,2,(kmer-1))]))/length(sequ[2])
        end
    end
    #Save the zval for this sequence
    zval[1:4^kmer,counter]=collect(values(sort(dictres)))

    #Increase the counter
    global counter+=1
end

#Compute the correlation matrix between z-val profiles
results=cor(zval,dims=1)
#results=corspearman(zval)

function wout(resmat,seqs)
    #Prepare results matrix for writing out
    resmat=string(resmat)
    ind1=findfirst("[",resmat)[1]+1
    ind2=findfirst("]",resmat)[1]-1
    resmat=resmat[ind1:ind2]
    resmat=replace(resmat,"; " => "\n")
    resmat=replace(resmat," " => "\t")

    #Write out results
    io=open("ResMat.txt","w")
    write(io,resmat)
    close(io)

    #Write out tags of the sequences
    tagio=open("SeqNames.txt","w")
    for l in seqs
        write(tagio,l[1])
        write(tagio,"\n")
    end
    close(tagio)
end

wout(results,sequences)
