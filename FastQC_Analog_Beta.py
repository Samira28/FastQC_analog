def gc_content(file, outdir):

    record_dict1 = SeqIO.to_dict(SeqIO.parse(file, "fastq"))
    GC = []
    for key in record_dict1.keys():
        G = []
        C = []
        N = []
        for letter in record_dict1[key]:
            if letter == "G":
                G.append("letter")
            if letter == "C":
                C.append("letter")
            if letter == "N":
                N.append("letter")
        a = ((len(G)+len(C))/(len(record_dict1[key]) - len(N)))*100
        GC.append(a)
    return sns.distplot(GC) 
    return sns_plot.figure.savefig(outdir + "output_gc_content.png")

def get_read_length(records):
    for record in records:
        return len(record.seq)

def error_probability(file, outdir):

    with open(file,'rt') as handle:
        num_reads = 0
        records = SeqIO.parse(handle, format = "fastq")
        len_read = get_read_length(records)
        p = [0]*len_read
        for i,record in enumerate(records):
            num_reads += 1
            quality = record.letter_annotations["phred_quality"]
            for pos, q in enumerate(quality):
                p[pos] += 10**(-quality[pos]/10)/(1+10**(-quality[pos]/10))
    for pos in range(len_read):
        p[pos] = float(p[pos])/num_reads

    plt.xlabel("Position")
    plt.ylabel("Error probability")
    plt.plot(range(len_read), p, color = "blue")
    return plt.show()
    return plt.savefig(outdir + 'error_probability.png')


def qual_distributon(file, outdir):

    with open(file,'rt') as handle:
        num_reads = 0
        records = SeqIO.parse(handle, format = "fastq")
        len_read = get_read_length(records)
        q = [0]*len_read
        for i,record in enumerate(records):
            num_reads += 1
            quality = record.letter_annotations["phred_quality"]
            for pos, t in enumerate(quality):
                q[pos] += quality[pos]
    for pos in range(len_read):
        q[pos] = float(q[pos])/num_reads

    plt.xlabel("Position")
    plt.ylabel("Quality")
    plt.plot(range(len_read), q, color = "blue")
    return plt.show()
    return plt.savefig(outdir + 'qual_distributon.png')

def distribution_insert_distance(file, outdir):

    samfile = pysam.AlignmentFile( "/Users/samira/IB_PRAC/Project_1/alignment.sam", "rb")
    tlen_and_suffering = []
    for read in samfile.fetch():
        if read.template_length < 1500 and read.template_length > 0:
            tlen_and_suffering.append(read.template_length)

    num_bins = 100
    plt.hist(tlen_and_suffering, num_bins,color = "red")
    plt.xlabel("The distance of insertion")
    plt.ylabel("Quantity of reads")
    plt.show()
    plt.savefig(outdir + 'distribution_insert_distance.png')


def fastqc_analog(request, input_file, outdir):
    while request != "exit":   
        if request == "GC content":      
            return gc_content(input_file, outdir)
            print("Enter command")
        elif request == "quality distribution":
            return qual_distributon(input_file, outdir)
            print("Enter command")
        elif request == "error probability":
            return error_probability(input_file, outdir)
            print("Enter command")
        elif request == "insert distribution":
            return distribution_insert_distance(input_file, outdir)
            print("Enter command")
        else:
            print("Use these comands: GC content, quality distribution, error probability, insert distribution")
        request = input()
    print("Good luck!")

print("Wellcome to FactQC Analog! There are only several of functions available, as it is a beta version. Have a good time using FastQC Analog beta and wait for final version:) Please, Enter the command:")
request = input()
input_file = "/Users/samira/IB_PRAC/Project_1/amp_res_1.fastq"
outdir = "/Users/samira/IB_PRAC/Project_1/FASTQC_RESULTS/"
fastqc_analog(request, input_file, outdir)
