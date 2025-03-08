import argparse
import datetime
import glob
import gzip
import os
import fnmatch
import random
import shutil
import string

from Bio import SeqIO
from Bio import AlignIO
from Bio import Seq

from itertools import islice
from Bio.Align import AlignInfo
from Bio.Blast.Applications import NcbiblastnCommandline
import pandas as pd

import matplotlib.pyplot as plt
from statistics import mean
#.............................................................................................................................
#TODO: Convert print strings/variables to fstrings
def Concat_Sequential_Nanopore_Files(Directory, OutPutDIR):
    """
    This function will create a new directory and store sequentially combined fastq.gz reads (e.g. file 1, file 1 + file 2, file 1 + file 2 + file 3)
    Params:
        Directory - str, Path to the fastq.gz files
        ProjectName - str, ID Unique to this project
        OutPutDIR - str, Path to the newly created directory that contains the sequential combined fastq.gz files
    Returns:
        NewDirectory - str, Path to the newly created directory that contains the sequential combined fastq.gz files
        CombinedFilesListUpdated - lst, list of sequentially combined fastq.gz files
    """
    # Create OutputDIR in Directory
    NewDirectory = Make_Directory(ProjectName = Directory, Suffix = OutPutDIR)

    # Count number of Nanopore Reads in Directory
    CountofFastqGZ = len(fnmatch.filter(os.listdir(Directory), '*.fastq.gz'))

    # Combine Reads if there are more than one Nanopore Read
    File_and_Path_List = glob.glob(Directory + '*.fastq.gz')
    File_and_Path_List = sorted(File_and_Path_List)
    
    if CountofFastqGZ > 1:
        # i is the indexing value, this loop will concat files 0 to i
        i = 1
        while i < len(File_and_Path_List) + 1:
            Combined_File = f'{i:03}' + "_Combined_Minion.fastq.gz"
            # on the first pass of the loop,  we want to copy the file at index 0
            if i == 1:
                print("Copy: " + File_and_Path_List[0] + " to " + NewDirectory, flush = True)
                Copy_File2(CurrentDIR="", NewDIR=NewDirectory, SingleFile=File_and_Path_List[0], NewFileName=Combined_File)
                print("Copy Complete: " + File_and_Path_List[0] + " to " + NewDirectory, flush = True)
            else:
                # Join the first i number of files
                Partial_List = File_and_Path_List[:i]
                CommandString = 'cat ' + ' '.join(Partial_List) + ' > ' + NewDirectory + Combined_File
                print(CommandString, flush = True)
                os.system(CommandString)
                print(str(i) + " complete", flush = True)
            i += 1
    else:
        print("Single File Detected", flush = True)
        print("Copy: " + File_and_Path_List[0] + " to " + NewDirectory, flush = True)
        Copy_File2(CurrentDIR="", NewDIR=NewDirectory, SingleFile=File_and_Path_List[0], NewFileName=Combined_File)
        print("Copy Complete: " + File_and_Path_List[0] + " to " + NewDirectory, flush = True)
    CombinedFileList = glob.glob(NewDirectory + '*.fastq.gz')
    CombinedFileListUpdated = [sub.replace('//', '/') for sub in CombinedFileList]
    CombinedFileListUpdated = sorted(CombinedFileListUpdated)
    return NewDirectory, CombinedFileListUpdated

def Convert_Fastq_Fasta(ProjectName, FastQ_Input, Subset, Suffix):
    """
    This function will uncompress and extracts reads from a single fastq.gz to a single fasta file. The number of reads that 
    are extracted from the fastq.gz file is set by the Subset variable.
    Params:
        ProjectName - str, ID Unique to this project
        FastQ_Input - str, fastq.gz file
        Subset - int, number of individual reads to extract
        Suffix - str, string that is added to the file name
    Returns:
        Fasta_File - str, fasta file that contains the extracted reads
    """
    Fasta_File = FastQ_Input.split(".")[0] + "_" + Suffix + ".fa"
    Fasta_File = ProjectName + "_" + Fasta_File.split("/")[-1]
    # make fastq
    if Subset != False:
        with gzip.open(FastQ_Input, "rt") as fastQ, open(Fasta_File, "w") as fastA:
            for record in islice(SeqIO.parse(fastQ, "fastq"), Subset):
                fastA.write(">" + record.id + "\n")
                fastA.write(str(record.seq) + "\n")
        fastA.close()
        fastQ.close()
    else:
        with gzip.open(FastQ_Input, "rt") as fastQ, open(Fasta_File, "w") as fastA:
            for record in SeqIO.parse(fastQ, "fastq"):
                #print(record.seq)
                fastA.write(">" + record.id + "\n")
                fastA.write(str(record.seq) + "\n")
        fastA.close()
        fastQ.close()
    print("Fasta file Created", flush = True)
    return Fasta_File

def Alignment_Muscle(Fasta_File, ProjectName):
    """
    This function calls the program muscle to perform an sequence alignment.
    Params:
        Fasta_File - str, fasta file
        ProjectName - str, ID Unique to this project
    Returns:
        Alignment_File - str, aligned sequences in a fasta file (.fa) format
    """
    Alignment_File = ProjectName + "_" + Fasta_File.split(".fa")[0] + "_subsetALN.fa"
    os.system("muscle -in " + Fasta_File + " -out " + Alignment_File)
    print("Alignment Complete", flush = True)
    return Alignment_File

def Extract_Alignment_Consensus(Alignment_File, Threshold, Require_Multiple, Header, ProjectName):
    """
    This function utilizes Bio.AlignIO from the BioPython library to extract a quick/rough consensus sequence from a sequence alignment
    Params:
        Alignment_File - str, fasta file
        Threshold - int, If the percentage of the most common residue type exceeds the set threshold, that residue type will output; 
            if not ambiguous character
        Header - str, header that will be written to the consensus fasta file
        ProjectName - str, ID Unique to this project
    Returns:
        consensus_fasta_file - str, consensus sequence of the alignment file in fasta (.fa) format
    """
    consensus_fasta_file = ProjectName + "_" + Alignment_File.split("_subsetALN.fa")[0] + "_subsetConsensus.fa"
    Alignment  = AlignIO.read(open(Alignment_File), "fasta")
    Consensus = AlignInfo.SummaryInfo(Alignment).gap_consensus(threshold=Threshold, ambiguous='N', require_multiple=Require_Multiple)
    Consensus = str(Consensus).replace('-','')
    with open(consensus_fasta_file, "w") as fastA:
        fastA.write(">" + Header + "\n" + Consensus + "\n")
    fastA.close
    print("Consensus Sequence complete", flush = True)
    return consensus_fasta_file

def PredictReference2(BlastInput, ProjectName, Suffix, BlastReference):
    """
    This function calls the program blastn to compare a sequence to reference fasta files.
    If blast fails because the rough consensus sequence generated in Extract_Alignment_Consensus was garbage,
    ASFV-Georgia will be utilized as an initial reference for downstream mapping
    Params:
        BlastInput - str, sequence in a fasta file to be analyzed
        ProjectName - str, ID Unique to this project
        Suffix - str, string to be added to the out file
        BlastReference - str, fasta file containing curated ASFV p72 sequences
    Returns:
        PredictedReference - str, The reference sequence that exhibits the highest percentage identity with the sampled sequence.
        BlastOutput - str, csv file containing results
    """
    BlastOutput = ProjectName + "_" + Suffix + "_blastoutput.csv"
    os.system(str(NcbiblastnCommandline(cmd='blastn', query = BlastInput, subject= BlastReference, max_hsps = 1, out = BlastOutput, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send", evalue=0.001)))
    # if the initial alignment fails, default to returning Georgia as Consensus so the script does not fail
    if Empty_File_Check(BlastOutput) == False:
        PredictedReference = pd.read_csv(BlastOutput, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send')).groupby('sseqid')['bitscore'].sum().sort_values(ascending=False).reset_index().iloc[0]['sseqid'] 
    else:
        PredictedReference = "ASFV_Georgia_2007_1"
    print("Blast Prediction Complete:" + PredictedReference, flush = True)
    return PredictedReference, BlastOutput

def Generate_Mapping_Stats_Simple(BAM_File, ProjectName, Suffix):
    """
    This function calls the utility mpileup from samtools to extract mapping quality and coverage at each position of the reference.
    *** 
    !!! The BAM file should have undergone prior processing, as seen in the MiniMapToReferenceSimple function.
    !!! It is essential that the BAM file exclusively includes sorted, mapped reads.
    ***
    Params:
        BAM_File - str, read alignment file in BAM format that had been previously massaged
        ProjectName - str, ID Unique to this project
        Suffix - str, string to be added to the out file
    Params (samtools mpileup):
        -s Output mapping qualities encoded as ASCII characters.
        -a Output all positions, including those with zero depth.
        -d At a position, read maximally INT reads per input file. Setting this limit reduces the amount of memory 
        and time needed to process regions with very high coverage. Passing zero for this option sets it to the 
        highest possible value, effectively removing the depth limit. [8000] Note that up to release 1.8, samtools 
        would enforce a minimum value for this option. This no longer happens and the limit is set exactly as specified.
    Returns:
        Final_Genome_Stats_File - str, tab separated text file. 'Genome', 'Position', 'Coverage', 'Quality_Cigar'

        ***
        TODO Note: this stat file is later edited in  the function Create_Graph2 where the mapping quality ASCII is converted to
        an average integer. This should probably be a function separate than Create_Graph2.
        ***
    """
    Final_Genome_Stats_File = ProjectName + "_" + Suffix + "_statistics.txt"
    os.system("samtools mpileup " + BAM_File + " -s -a -d 0 | cut -f3,5,6 --complement > " + Final_Genome_Stats_File)
    print("***_Mapping Stats 2/2 Complete_***", flush = True)
    return Final_Genome_Stats_File

def Rename_Header(fasta_file, string):
    """
    This function utilizes SeqIO from the BioPython library the program muscle to rename the header of
    a fasta file. This function should only be used with fasta files containing a single header/sequence.
    Params:
        fasta_file - str, fasta file with header to be renamed
        string - str, new header
    """
    with open(fasta_file, "r") as fastA:
        for record in SeqIO.parse(fastA, "fasta"):
            #new_header_and_sequence = ">" + string + "\n" + str(record.seq) + "\n"
            new_header_and_sequence = f">{string}\n{str(record.seq)}\n"
    fastA.close()
    with open(fasta_file, "w") as fastA:
        fastA.write(new_header_and_sequence)
    fastA.close()

def Empty_File_Check(file):
    """
    This function checks to see if a file is empty, does not exist, or if the variable
    is empty. The appropriate message is printed. Debugger tool
    Param:
        file - str, file to be examined
    Returns:
        messaged printed to console
    """
    if file != None:
        if os.path.exists(file):
            File_Size = os.stat(file).st_size
            if(File_Size == 0):
                print(f"{file} has zero data", flush = True)
                return True
            else:
                return False
        else:
            print(f"{file} does not exist", flush = True)
            return True
    else:
        print(f"{file} is an empty variable", flush = True)
        return True

def Combine_Stats(ProjectName, List_of_stats):
    """
    This function concats multiple tab separated stats files into a two pandas dataframes. One data
    frame contains all coverage stats at each position of the genome while the other dataframe is
    composed of all quality mapping stats at each position of the genome. 
    Params:
        ProjectName - str, ID Unique to this project
        List_of_stats - lst, list of stat files (tab separated text files) generated in "Generate_Mapping_Stats_Simple"
        AND converted in Create_Graphs2
    Returns:
        Coverage_Graph - str, graph displaying the coverage at each genome position (.png)  
        Quality_Graph - str, graph displaying the mapping quality at each genome position (.png)  
    """
    Final_Coverage_DF = pd.DataFrame(columns = ['Position', 'Coverage', 'Mapping_Quality']).set_index('Position')
    Final_Mapping_DF = pd.DataFrame(columns = ['Position', 'Mapping_Quality']).set_index('Position')
    Column_List = []
    for file in List_of_stats:
        SampleID = file.split("_")[4]
        Column_List.append(SampleID)
        Current_DF = pd.read_csv(file, sep = "\t", quoting=3).set_index('Position')
        Current_DF = Current_DF.drop(columns=['Genome', 'Quality_Cigar'])
        Final_Coverage_DF[SampleID] = Current_DF['Coverage']
        Final_Mapping_DF[SampleID] = Current_DF['Mapping_Quality']
    
    # Graph Coverage and Save PNG
    plt.close("all")
    CoverageGraph = Final_Coverage_DF.plot(kind = 'line',
        use_index=True,
        y = Column_List,
        logy = True,
        legend = True,
        title = 'Read Coverage',
        ylabel = 'Depth of Coverage',
        figsize=(20, 10)
    )
    CoverageFigure = CoverageGraph.get_figure()
    CoverageFigure.savefig(ProjectName + '_Compiled_Graph.png')
    plt.close("all")
    # Graph Quality and Save PNG
    QualityGraph = Final_Mapping_DF.plot(kind = 'line',
        use_index=True,
        y = Column_List,
        logy = False,
        legend = True,
        title = 'Mapping Quality',
        ylabel = 'Mapping Quality (max 60)',
        figsize=(20, 10)
    )
    CoverageFigure = QualityGraph.get_figure()
    CoverageFigure.savefig(ProjectName + '_Compiled_QualityGraph.png')
    plt.close("all")
    return ProjectName + '_Compiled_Graph.png', ProjectName + '_Compiled_QualityGraph.png'

def Combine_Blastoutputs(ProjectName, List_of_Blast_Files, Blast_Type):
    """
    This function concats the top hit of all created blast files into a single file using pandas
    Params:
        ProjectName - str, ID Unique to this project
        List_of_Blast_Files - lst, list of blast csv files
        Blast_Type - str, string that is added to the out file
    Returns:
        FileOutput - str, file (csv) containing the top hits of all concatenated blast files
    """
    # set header names
    colnames =['qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send']
    # create empty dataframe
    Final_Blast_DF = pd.DataFrame(columns = colnames)
    # import first row from all blastoutputs, append row 1 to final DF, top hit should be row 1
    for file in List_of_Blast_Files:
        Current_DF = pd.read_csv(file, sep = "\t", quoting=3, header = None)
        Current_DF.columns = colnames
        First_Row = Current_DF.iloc[0]
        #print(First_Row)
        #Final_Blast_DF = pd.concat([Final_Blast_DF, First_Row], axis=1, ignore_index= True)
        Final_Blast_DF.loc[len(Final_Blast_DF.index)] = First_Row
    # create final file
    FileOutput = ProjectName + "_consensus_sequence_" + Blast_Type + "_summary.csv"
    Final_Blast_DF.to_csv(FileOutput, sep = ",", index = False)
    return FileOutput

def Combine_Consensus_Fastas(ProjectName, List_of_fastas):
    """
    This function concats multiple fasta files into a single fasta (.fa) file. 
    Each sequence retains its header
    Params:
        ProjectName - str, ID Unique to this project
        List_of_fastas - lst, list of fasta files
    Returns:
        Final_File - str, combined fasta file (.fa) 
    """
    string = ' '.join(List_of_fastas)
    Final_File = ProjectName + "_consensus_all.fa"
    os.system("cat " + string + " > " + Final_File)
    return Final_File

def MetaDataUpdate_batch(MetaDataFile, NewMetaData, ProjectName, Row):
    """
    This function updates the MetaDataStorage.csv file with the current run(s) info
    Params:
        MetaDataFile - str, MetaDataStorage.csv
        NewMetaData - str, metadata file (.csv) containing info on the current project and run
        ProjectName - str, ID Unique to this project
        Row - int,  run # based on the NewMetaData file. Determined from the function Number_of_Rows_MetaDataFile
    Returns:
        updates MetaDataStorage.csv
    """
    MetaDataDF = pd.read_csv(MetaDataFile)
    MetaDataDF = MetaDataDF.where(pd.notnull(MetaDataDF), None)
    
    NewMetaDF = pd.read_csv(NewMetaData).iloc[Row]
    NewMetaDF["Project_ID (internal)"] = ProjectName

    MetaDataToStore = pd.concat([MetaDataDF,NewMetaDF.to_frame().transpose()],axis=0)
    NewMetaDF = NewMetaDF.where(pd.notnull(NewMetaDF), None)

    MetaDataToStore.to_csv(MetaDataFile, index=False)

    print("***_Metadata Updated_***", flush = True)
    return NewMetaDF

def Number_of_Rows_MetaDataFile(NewMetaData):
    """
    When running this pipeline as the batch, this function will determine the number of samples/ runs
    Params:
        NewMetaData - str, csv file containing the sample/runs metadata
    Returns:
        Number_of_rows - int, number of samples to be analyzed when the pipeline is initiated
    """
    NewMetaDF = pd.read_csv(NewMetaData)
    Number_of_rows = len(NewMetaDF)
    return Number_of_rows

def Consensus_Extractor_p72pipeline(Genome_Fasta, BAM_File, Max_Depth_Coverage, Suffix, ProjectName):
    """
    This function conducts an analysis of a read mapping file (BAM) in relation to the indexed reference 
    genome to which the reads were initially aligned. Utilizing the mpileup utility from bcftools, it 
    generates a multi-way pileup of genotype likelihoods, which is subsequently piped into the bcftools 
    call to produce a VCF file that encompasses all identified variants. This VCF file is then transformed 
    into a pandas dataframe, where each variant is analyzed to if it appears in 50% or more of the reads
    at a specific position. Variants that meet this criterion are preserved. Ultimately, a final VCF is 
    generated, which is then compressed and indexed using bcftools view.
    *** 
    !!! The BAM file should have undergone prior processing, as seen in the MiniMapToReferenceSimple function.
    !!! It is essential that the BAM file exclusively includes sorted, mapped reads.
    ***
    Params:
        Genome_Fasta - str, Reference genome fasta file
        BAM_File - str, read mapping file (BAM) previously massaged
        Max_Depth_Coverage, str - see bcftools mpileup -L
        Suffix - str, string to be added to the out file
        ProjectName - str, ID Unique to this project
    Params (bcftools mpileup)
        -L - int, Skip INDEL calling if the average per-sample depth is above INT [250]
        -Ou -str, output type, Use the -Ou option when piping between bcftools subcommands 
            to speed up performance by removing unnecessary compression/decompression 
            and VCF←→BCF conversion
        -f - str, The faidx-indexed reference file in the FASTA format. The file can be 
            optionally compressed by bgzip. Reference is required by default unless the 
            --no-reference option is set [null]
        -d - int, At a position, read maximally INT reads per input file. Note that the 
            original samtools mpileup command had a minimum value of 8000/n where n was the 
            number of input files given to mpileup. This means that in samtools mpileup the 
            default was highly likely to be increased and the -d parameter would have an 
            effect only once above the cross-sample minimum of 8000. This behavior was 
            problematic when working with a combination of single- and multi-sample bams, 
            therefore in bcftools mpileup the user is given the full control 
            (and responsibility), and an informative message is printed instead [250]
    Params (bcftools call)
        -mv -m alternative model for multiallelic and rare-variant calling designed to overcome 
            known limitations in -c calling model (conflicts with -c)
            -v output variant sites only
        -Ov -str, output type -Ov option for uncompressed VCF
        -o- str, output
    Params (bcftools view)
        -Oz compressed VCF 
        --write-index Automatically index the output file
        -o out file
    Returns:
        Variant_Table - str, User friendly file (csv) which displays variant info
        VCF_Output
        VCF_index

    """
    #Find all SNPs
    #
    RAW_File = ProjectName + "_" + Suffix + "_calls.vcf"
    #os.system("bcftools mpileup -Ou -f " + Genome_Fasta + " " + BAM_File + " -d " + Max_Depth_Coverage + " | bcftools call -mv -Ov -o "  + RAW_File)
    os.system("bcftools mpileup -L " + Max_Depth_Coverage + " -Ou -f " + Genome_Fasta + " " + BAM_File + " -d " + Max_Depth_Coverage + " | bcftools call -mv -Ov -o "  + RAW_File)
    print("***_Consensus Extractor 1/3 Complete_***", flush = True)
    ###############################################################################
    #
    # the first 27 rows can be skipped because they are not formatted properly
    Variant_Table = pd.read_csv(RAW_File, sep="\t", skiprows = 27)
    RAW_Table = pd.read_csv(RAW_File, sep="\t", skiprows = 27)
    
    # Stop Process if there are no Variants, 
    # still need to run bcftools to create index file to correctly run apply_variant function
    if len(Variant_Table.index) == 0:
        VCF_index = RAW_File + ".gz"
        CSV_Output = ProjectName + "_" + Suffix + "_SNP.csv"
        print("no variants found", flush = True)
        os.system("bcftools view "+ RAW_File + " -Oz --write-index -o " + VCF_index)
        return VCF_index, CSV_Output
    """if len(Variant_Table.index) == 0:
        VCF_index = ''
        print("no variants found")
        return VCF_index"""
    
    # set file names
    CSV_Output = ProjectName + "_" + Suffix + "_SNP.csv"
    VCF_Output = ProjectName + "_" + Suffix + "_corrected_calls.vcf"
    VCF_index = VCF_Output + ".gz"
    
    # the INFO column has an incredible amount of data that is separated by ';'
    Variant_Table['INFO'] = Variant_Table['INFO'].str.split(';').fillna(Variant_Table['INFO'])
    # explode the data into multiple rows but keep the original index # on each created row
    Variant_Table=Variant_Table.explode('INFO',ignore_index=False)

    # of the data in info, we only care about 'DP' (total coverage) and DP4 (Ref Read F, Ref Read R, Variant Read F, Variant Read R)
    searchfor = ['DP']
    Variant_Table = Variant_Table[Variant_Table.INFO.str.contains('|'.join(searchfor))]

    # Create another dataframe where DP and DP4 is separated using '=', name the columns ID and Value
    SplitColumn = Variant_Table['INFO'].str.split(pat = '=', expand = True)
    SplitColumn.columns = ['ID', 'Value']

    # Create DP and DP4 columns, since indices were preserved, the values for DP and DP4 at the same position will be preserved
    SplitColumn = SplitColumn.pivot(columns = 'ID', values = 'Value')

    # Restore string to integer
    CoverageColumn = SplitColumn['DP'].astype(int)

    # since DP4 is composed of (Ref Read F, Ref Read R, Variant Read F, Variant Read R), we only care about Variant Read F & R
    DP4Column = SplitColumn['DP4'].str.split(',', expand = True).astype(int)

    # The sum of this value gives us the number of reads that contain the variant
    VariantSum = DP4Column[2] + DP4Column[3]

    # The frequency of reads that compare the variants as compared to the total coverage at that position
    Frequency = round(VariantSum / CoverageColumn, 3)
    
    # Clean Up Table
    DropColumns = ['ID', 'FILTER', 'INFO', 'FORMAT', BAM_File]
    Variant_Table = Variant_Table[~Variant_Table.index.duplicated(keep='first')].drop(DropColumns, axis = 1)
    Variant_Table.rename(columns={
        "#CHROM": "Reference_Strain",
        "POS": "Position",
        "REF": "Reference",
        "ALT": "Sequenced_Sample"
        }
    )
    # Create columns with calculated values
    Variant_Table['Count'] = VariantSum
    Variant_Table['Coverage'] = CoverageColumn
    Variant_Table['Frequency'] = Frequency
    
    # Drop Variants less than 50%
    Variant_Table = Variant_Table[Variant_Table.Frequency >= 0.50]
    
    # Drop indices in original file that are not in the final SNP Table
    RAW_Table = RAW_Table[RAW_Table.index.isin(Variant_Table.index)]

    # Save Variant File in an User Friendly Format
    Variant_Table.to_csv(CSV_Output, sep=',', index=False)
    
    # Creating A Corrected Variant File for PipeLine
    with open(RAW_File) as myfile:
        first_27_lines = myfile.readlines()[0:27]
    myfile.close()
    with open(VCF_Output,'a') as Corrected_File:
        for line in first_27_lines:
            Corrected_File.write(line)
    Corrected_File.close()
    RAW_Table.to_csv(VCF_Output, sep = '\t', mode='a', index = False, header=True)
    #
    print("***_Consensus Extractor 2/3 Complete_***", flush = True)
    ##############################################################################
    os.system("bcftools view "+ VCF_Output + " -Oz --write-index -o " + VCF_index)
    print("***_Consensus Extractor 3/3 Complete_***", flush = True)
    return VCF_index, CSV_Output

def Make_Directory(ProjectName, Suffix):
    """
    This function creates a new directory using the ProjectName and Suffix stings if it does not exist.
    Params:
        ProjectName - str, ID Unique to this project
        Suffix - str, string to be added to the out file
    Returns:
        NewDirectory - str, path to the new directory
    """
    NewDirectory = ProjectName + Suffix + "/"
    if not os.path.exists(NewDirectory):
        os.makedirs(NewDirectory)
        return NewDirectory
    else:
        return NewDirectory

def Copy_File2(CurrentDIR, NewDIR, SingleFile, NewFileName):
    """
    This function copies/renames a file to a new directory.
    Params:
        CurrentDIR - str, ID Unique to this project
        NewDIR - str, string to be added to the out file
        SingleFile - str,
        NewFileName - str,
    """
    if not os.path.exists(NewDIR):
        os.makedirs(NewDIR)
    else:
        pass

    if os.path.exists(CurrentDIR + SingleFile):
        shutil.copy2(CurrentDIR + SingleFile, NewDIR + NewFileName)
    else:
        print("Error, no file to copy", flush = True)

def MiniMapToReferenceSimple(GenomeFileTemp, CombinedMinionReadFile, ProjectName, Suffix):
    """
    This function uses minimap2 to align nanopore reads to a reference genome, 
    generating a SAM file. The SAM file is then filtered and converted to a BAM file 
    with samtools view, removing unmapped reads. Finally, the reads are sorted by 
    genomic position using samtools sort.
    Params:
        GenomeFileTemp - str, fasta file to which the reads will be aligned to
        CombinedMinionReadFile - Read fastq.gz file
        ProjectName - str, ID Unique to this project
        Suffix - str, string to be added to the out files
    Params (minimap2):
        -a Generate CIGAR and output alignments in the SAM format.
        -x map-ont preset option. Align noisy long reads of ~10% error rate to a reference genome.
    Params (samtools view):
        -b Output in the BAM format.
        -F 4 Do not output alignments with any bits set in FLAG present in the FLAG field
            The FLAG preset 4 indicates "read unmapped" -> "Do not output unmapped reads"
    Returns:
        list [
            MinionMappedTempSAM - temp .sam file,
            MinionMappedTempBAM - temp .bam file,
            MinionMappedTempSOR - read alignment file (.bam)
        ]
    """
    MinionMappedTempSAM = ProjectName + "_" + Suffix + "_Minion_Mapped.sam"
    MinionMappedTempBAM = ProjectName + "_" + Suffix + "_Minion_Mapped.bam"
    MinionMappedTempSOR = ProjectName + "_" + Suffix + "___Minion_Sort.bam"
    #os.system("minimap2 -ax map-ont " + GenomeFileTemp + " " + CombinedMinionReadFile + " -o " + MinionMappedTempSAM)
    os.system("minimap2 -ax map-ont " + GenomeFileTemp + " " + CombinedMinionReadFile + " > " + MinionMappedTempSAM)
    print("***_MiniMap 1/3 Complete_***", flush = True)
    os.system("samtools view -b -F 4 " + MinionMappedTempSAM + " > " + MinionMappedTempBAM)
    print("***_MiniMap 2/3 Complete_***", flush = True)
    os.system("samtools sort " + MinionMappedTempBAM + " -o " + MinionMappedTempSOR)
    print("***_MiniMap 3/3 Complete_***", flush = True)
    #os.system("samtools index " + MinionMappedTempSOR)
    #print("***_MiniMap 4/6 Complete_***")
    return [MinionMappedTempSAM, MinionMappedTempBAM, MinionMappedTempSOR]

def Delete_File(DeleteFile):
    """
    This function will delete a file if it exists. . An error message is returned if
    the file does not exist.
    Params:
        DeleteFile - str, file to be deleted
    """
    if DeleteFile != None:
        if os.path.exists(DeleteFile):
            os.remove(DeleteFile)
        else:
            print(str(DeleteFile + "does not exist"), flush = True)
    else:
        print("empty variable", flush = True) 

def UnMappedRegions(BAM_File, Suffix, ProjectName):
    """
    This function utilizes bedtools genomecov to  determine if there are 
    regions on the reference genome where no reads align.
    *** 
    !!! The BAM file should have undergone prior processing, as seen in the MiniMapToReferenceSimple function.
    !!! It is essential that the BAM file exclusively includes sorted, mapped reads.
    ***
    Params:
        BAM_File - str, read-alignment file
        Suffix - str, string to be added to the out file
        ProjectName - str, ID Unique to this project
    Params (bedtools genomecov)
        -ibam	BAM file as input for coverage. Each BAM alignment in A added to the total coverage for the genome.
        -bga	Report depth in BedGraph format. Regions with zero coverage are also reported. This allows one to 
                quickly extract all regions of a genome with 0 coverage by applying: “grep -w 0$” to the output.
    """
    UnMappedRegions_Bed_File = ProjectName + "_unmappedRegions_" + Suffix + ".bed"
    os.system("bedtools genomecov -ibam " + BAM_File + " -bga | grep -w 0$ > " + UnMappedRegions_Bed_File)
    print("***_UnMapped Regions Complete_***", flush = True)
    return UnMappedRegions_Bed_File

def Apply_Variants_to_Consensus(Genome_Fasta, VCF_Index, BED_File, ProjectName, Suffix):
    """
    This function uses bcftools consensus to incorporate variants and unmapped regions into 
    a consensus genome. Essentially, once the reads have been aligned to a reference genome 
    (MiniMapToReferenceSimple) and variants that are supported by more than 50% of the reads 
    are identified (Consensus_Extractor_p72pipeline), these variants are then integrated into 
    the reference genome. Furthermore, any regions that lack coverage (UnMappedRegions) are 
    designated as 'N' to indicate ambiguous nucleotides.
    Params:
        Genome_Fasta - str, reference genome
        VCF_Index - str, indexed VCF variant file (see function Consensus_Extractor_p72pipeline)
        BED_File - str, bed file (see function UnMappedRegions)
        ProjectName - str, ID Unique to this project
        Suffix - str, string to be added to the out files
    Params (bedtools genomecov)
        -ibam	BAM file as input for coverage. Each BAM alignment in A added to the total coverage for the genome.
        -bga	Report depth in BedGraph format. Regions with zero coverage are also reported. This allows one to 
                quickly extract all regions of a genome with 0 coverage by applying: “grep -w 0$” to the output.
    Returns
        Consensus_Sequence - str, fasta file (.fa) of the read-aligned consensus sequence
    """
    Consensus_Sequence = ProjectName + "_" + Suffix + ".fa"
    os.system("cat " + Genome_Fasta + " | bcftools consensus " + VCF_Index + " -m " + BED_File + " > " + Consensus_Sequence)
    print("***_Apply Variants to Fasta Complete_***", flush = True)
    return Consensus_Sequence

def Create_Graph2(Stats_File, ProjectName, Suffix):
    """
    This function uses the pandas and matlab library and is designed to generate a coverage 
    graph alongside a graph that illustrates the quality of read mapping, utilizing a 
    previously produced tab-separated statistics text file 
    (refer to the function Generate_Mapping_Stats_Simple). This file comprises the 
    headers: Genome, Position, Coverage, and Quality Cigar. Additionally, the 
    function processes the ASCII Phred 33 characters found in the Quality_Cigar 
    column, converting and averaging them into integer values. These computed values 
    are then recorded in a newly established Mapping_Quality column, and the 
    statistics file is subsequently updated and saved to incorporate this additional column.
    TODO: Create a separate function to convert ascii
    Params:
        Stats_File - str, tab separated text file (see function Generate_Mapping_Stats_Simple)
        ProjectName - str, ID Unique to this project
        Suffix - str, string to be added to the out files
    Returns:
        CoverageGraph - str, graph file (.png) demonstrating read coverage across reference genome
        QualityGraph - str, , graph file (.png) demonstrating read mapping qc across reference genome
        AverageCoverage - int, average read coverage across the genome
    """
    plt.close("all")
    StatDataFrame = pd.read_csv(Stats_File, sep='\t', header = None, names = ['Genome', 'Position', 'Coverage', 'Quality_Cigar'], quoting=3)
    
    # Graph Coverage and Save PNG
    CoverageGraph = StatDataFrame.plot(kind = 'line',
        x = 'Position',
        y = 'Coverage',
        logy = True,
        legend = False,
        title = 'Read Coverage',
        ylabel = 'Depth of Coverage',
        figsize=(20, 10)
    )
    CoverageFigure = CoverageGraph.get_figure()
    CoverageFigure.savefig(ProjectName + "_" + Suffix + '_CoverageGraph.png')
    plt.close("all")
    print("CoverageFigure Done", flush = True)

    # Calculate Average Coverage
    AverageCoverage = round(StatDataFrame.Coverage.mean(), 2)

    # Graph Mapping Quality and Save PNG
    StatDataFrame.insert(3, 'Mapping_Quality', StatDataFrame['Quality_Cigar'].apply(lambda x: round(mean([ord(i) for i in x]) - 33, 2)))
    

    
    QualityGraph = StatDataFrame.plot(kind = 'line',
        x = 'Position',
        y = 'Mapping_Quality',
        legend = False,
        title = 'Read Mapping Quality',
        ylabel = 'Mapping Quality (max 60)',
        figsize=(20, 10)
    )
    QualityFigure = QualityGraph.get_figure()
    QualityFigure.savefig(ProjectName + "_" + Suffix  + '_QualityGraph.png')

    # Update CSV
    StatDataFrame.to_csv(Stats_File, sep='\t', index = False)


    print("***_Graphs Complete_***", flush = True)
    return ProjectName + "_" + Suffix  + '_CoverageGraph.png', ProjectName + "_" + Suffix   + '_QualityGraph.png', AverageCoverage

def rand_pass(size):
    """
    This function generates a random string of characters and digits. This will be used
    as a unique ID to create the ProjectName.
    Param:
        size - int, length of the generated string
    Return:
        generate_pass - str, random string of n characters and digits
    """
    generate_pass = ''.join([random.choice(string.ascii_lowercase + string.digits) for n in range(size)]) 
    return generate_pass

def Move_File(Date, ProjectName, List_of_Files):
    """
    This function moves each file in a list to a new directory that is named for the current project.
    TODO: change Date variable to something more relevant (ie sample name)
    Params:
        Date - str, Prefix added to the Project directory
        ProjectName - str, ID Unique to this project
        List_of_Files - lst, list of file to be moved
    """
    NewDirectory = Date + "_Project_" + ProjectName +"/"
    if not os.path.exists(NewDirectory):
        os.makedirs(NewDirectory)
    else:
        pass
    for file in List_of_Files:
        if file != None:
            if os.path.exists(file):
                shutil.move(file, NewDirectory + file)
        else:
            print("empty variable", flush = True)

#.............................................................................................................................
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Setting Variables ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
#.............................................................................................................................

# Files Present in the current directory
P72_NUCLEOTIDE_DB = "p72_References/"
BLASTREFERENCE = "p72_References/054_DNA_ND_p72.fa"
METADATASTORAGE= "MetaDataStorage.csv"

#..........................................................................................
parser = argparse.ArgumentParser(description='Mapping Nanopore reads to ASFV P72')
parser.add_argument('--metadata', action="store", dest='metadata', default=0)
args = parser.parse_args()
NEWMETADATA = args.metadata
if NEWMETADATA == 0:
    print("please indicate the metadata file using  --metadata \n Example python3 P72_arg_Nanopore_Script --metadata MetadataNew.csv", flush=True)
    quit()
#..........................................................................................


# SubSampling Variables
MaxSubSample = 10
Threshold = 0.2
Require_Multiple = 1
Sample_Counter = 0

# Determine number of Samples
Number_of_Samples = Number_of_Rows_MetaDataFile(NEWMETADATA)


# This is the main Loop. It will run the process for every sample in "NEWMETADATA"...............................................
while Sample_Counter < Number_of_Samples:

    #.............................................................................................................................
    # ~ ~ ~ ~ ~ ~ ~ ~ Creating a quick consensus from Minion Reads and determining a rough reference  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    #.............................................................................................................................

    # Create Project ID
    TIME = str(datetime.date.today())
    ProjectName = TIME + "_" + rand_pass(8)
    ProjectName = ProjectName.replace(" ", "_").replace("-","_")

    # Update MetadataArchive, Set Nanopore Directory
    NewMetaDF = MetaDataUpdate_batch(MetaDataFile = METADATASTORAGE, NewMetaData = NEWMETADATA, ProjectName = ProjectName, Row = Sample_Counter)
    Directory = NewMetaDF['MinionDirectory']

    # Combine Nanopore Files Sequentially
    Combined_Minion_Directory, Combined_Minion_List = Concat_Sequential_Nanopore_Files(Directory=Directory, OutPutDIR = ProjectName)
    MinionFiletoSample = Combined_Minion_List[0]

    # Create a Fasta file from a subset of fastq files
    Minion_Fasta_Subset = Convert_Fastq_Fasta(ProjectName=ProjectName, FastQ_Input = MinionFiletoSample, Subset = MaxSubSample, Suffix="subset")

    # Create an Alignment from the fasta file
    Alignment_File = Alignment_Muscle(Fasta_File = Minion_Fasta_Subset, ProjectName =ProjectName)

    # Extract Consensus from Alignment
    Temp_subset_consensus_fasta_file = Extract_Alignment_Consensus(Alignment_File = Alignment_File, Threshold = Threshold, Require_Multiple = Require_Multiple, Header="Consensus", ProjectName =ProjectName)

    # Predict the best reference using blastn, remove temp files
    Predicted_Course_Reference, BlastOutput_tempConsensus_File = PredictReference2(BlastInput = Temp_subset_consensus_fasta_file, 
    ProjectName = ProjectName, 
    Suffix = "Extracted_subset", 
    BlastReference = BLASTREFERENCE)
    p72_Genome_Course_File = P72_NUCLEOTIDE_DB + Predicted_Course_Reference + ".fa"
    #Delete_File(BlastOutput_tempConsensus_File)
    #Delete_File(Alignment_File)
    #Delete_File(Temp_subset_consensus_fasta_file)
    #Delete_File(Minion_Fasta_Subset)

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Complete Rough Reference ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush = True)

    #.............................................................................................................................
    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Mapping Nanopore Reads to Best Reference to Determine a Better Reference ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    #.............................................................................................................................

    # Variables
    SampleID = "temp_refine"

    # Map Nanopore Reads to Reference
    MappingOutputList_refine = MiniMapToReferenceSimple(GenomeFileTemp = p72_Genome_Course_File, CombinedMinionReadFile = Combined_Minion_List[0], ProjectName = ProjectName, Suffix = SampleID)
    # OutputList = [MinionMappedTempSAM, MinionMappedTempBAM, MinionMappedTempSOR]

    # Determine SNP
    Variant_File_refine, SNP_csv_file = Consensus_Extractor_p72pipeline(Genome_Fasta = p72_Genome_Course_File, BAM_File = MappingOutputList_refine[2], Max_Depth_Coverage = "10000", Suffix = SampleID, ProjectName = ProjectName)

    # Determine UnMapped Regions
    UnMapped_Bed_File_refine = UnMappedRegions(BAM_File = MappingOutputList_refine[2], Suffix = SampleID, ProjectName = ProjectName)

    # Create Consensus
    Temp_Sample_Consensus_p72_File = Apply_Variants_to_Consensus(Genome_Fasta = p72_Genome_Course_File, VCF_Index = Variant_File_refine, BED_File = UnMapped_Bed_File_refine, ProjectName = ProjectName, Suffix = SampleID + "_consensus")

    # Rename Header in Consensus
    new_header = ProjectName + "_" + SampleID + "_consensus"
    Rename_Header(fasta_file = Temp_Sample_Consensus_p72_File, string = new_header)

    # Determine True Reference:
    Predicted_Refined_Reference, Sample_BlastOutput_File_refine = PredictReference2(BlastInput = Temp_Sample_Consensus_p72_File, 
    ProjectName = ProjectName, 
    Suffix=SampleID, 
    BlastReference = BLASTREFERENCE)
    Final_p72_Reference = P72_NUCLEOTIDE_DB + Predicted_Refined_Reference + ".fa"

    # Delete Temp Files
    Delete_File(ProjectName + "_" + SampleID + "_calls.vcf")
    Delete_File(ProjectName + "_" + SampleID + "_calls.vcf.gz")
    Delete_File(ProjectName + "_" + SampleID + "_calls.vcf.csi")
    Delete_File(ProjectName + "_temp_refine_calls.vcf.gz.csi")
    Delete_File(ProjectName + "_" + SampleID + "_corrected_calls.vcf")
    Delete_File(ProjectName + "_" + SampleID + "_corrected_calls.vcf.gz")
    Delete_File(ProjectName + "_" + SampleID + "_corrected_calls.vcf.gz.csi")
    Delete_File(MappingOutputList_refine[0])
    Delete_File(MappingOutputList_refine[1])
    Delete_File(MappingOutputList_refine[2])
    Delete_File(Sample_BlastOutput_File_refine)
    Delete_File(UnMapped_Bed_File_refine)
    Delete_File(Temp_Sample_Consensus_p72_File)
    Delete_File(SNP_csv_file)
    del(SampleID)

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Complete Refined Reference ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush = True)
    #.............................................................................................................................
    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Mapping Nanopore Reads to Best Reference to Determine a Better Reference ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    #.............................................................................................................................

    # Creating Empty lists to store file names 
    Consensus_Fasta_List = []
    BlastN_Output___List = []
    Statistics_File_List = []
    Qual_Graph_File_List = []
    Cove_Graph_File_List = []
    BAM________File_List = []
    SNP_csv____File_List = []

    # Only Look at first 50 files
    if len(Combined_Minion_List) > 50:
        Combined_Minion_List = Combined_Minion_List[:50]
    # Loop for each Concat Nanopore File
    for Minion_Combo_Files in Combined_Minion_List:
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Analyzing " + Minion_Combo_Files + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush = True)
        #Current Sample Number
        SampleID = Minion_Combo_Files.split("fastq.gz")[0]
        SampleID = SampleID.split("/")[-1]
        SampleID = SampleID.split("_")[0]

        # Map Nanopore Reads to Reference
        MappingOutputList = MiniMapToReferenceSimple(GenomeFileTemp = Final_p72_Reference, CombinedMinionReadFile = Minion_Combo_Files, ProjectName = ProjectName, Suffix = SampleID)
        # OutputList = [MinionMappedTempSAM, MinionMappedTempBAM, MinionMappedTempSOR]
        # Delete Temp Files
        Delete_File(MappingOutputList[0])
        Delete_File(MappingOutputList[1])

        # Determine SNP
        Variant_File, SNP_csv_file = Consensus_Extractor_p72pipeline(Genome_Fasta = Final_p72_Reference, BAM_File = MappingOutputList[2], Max_Depth_Coverage = "10000", Suffix = SampleID, ProjectName = ProjectName)

        # Determine UnMapped Regions
        UnMapped_Bed_File = UnMappedRegions(BAM_File = MappingOutputList[2], Suffix = SampleID, ProjectName = ProjectName)

        # Calculate Stats
        BestHit_Stats_File = Generate_Mapping_Stats_Simple(BAM_File = MappingOutputList[2], ProjectName = ProjectName, Suffix = SampleID)

        # Create QC and Coverage Graphs
        Coverage_Graph_File, Quality_Graph_File, AverageCoverage = Create_Graph2(Stats_File = BestHit_Stats_File, ProjectName = ProjectName, Suffix=SampleID)

        # Create Consensus
        Sample_Consensus_p72_File = Apply_Variants_to_Consensus(Genome_Fasta = Final_p72_Reference, VCF_Index = Variant_File, BED_File = UnMapped_Bed_File, ProjectName = ProjectName, Suffix = SampleID + "_consensus")

        # Rename Header in Consensus
        new_header = ProjectName + "_" + SampleID + "_consensus"
        Rename_Header(fasta_file = Sample_Consensus_p72_File, string = new_header)

        # Determine Historic Genotype
        PredictedReference, Sample_BlastOutput_File = PredictReference2(BlastInput = Sample_Consensus_p72_File, ProjectName = ProjectName, Suffix=SampleID, BlastReference = BLASTREFERENCE)

        # Append Files to List
        Consensus_Fasta_List.append(Sample_Consensus_p72_File)
        BlastN_Output___List.append(Sample_BlastOutput_File)
        Statistics_File_List.append(BestHit_Stats_File)
        Qual_Graph_File_List.append(Quality_Graph_File)
        Cove_Graph_File_List.append(Coverage_Graph_File)
        BAM________File_List.append(str(MappingOutputList[2]))
        SNP_csv____File_List.append(SNP_csv_file)

        # Delete Temp Files
        Delete_File(ProjectName + "_" + SampleID + "_calls.vcf")
        Delete_File(ProjectName + "_" + SampleID + "_calls.vcf.gz")
        Delete_File(ProjectName + "_" + SampleID + "_calls.vcf.gz.csi")
        Delete_File(ProjectName + "_" + SampleID + "_corrected_calls.vcf")
        Delete_File(ProjectName + "_" + SampleID + "_corrected_calls.vcf.gz")
        Delete_File(ProjectName + "_" + SampleID + "_corrected_calls.vcf.gz.csi")
        Delete_File(UnMapped_Bed_File)

        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Complete " + Minion_Combo_Files + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush = True)

    #.............................................................................................................................
    # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ Final File Generation ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    #.............................................................................................................................
    OriginalSampleName = NewMetaDF['isolate']
    

    Final_Graph_File_1, Final_Graph_File_2 = Combine_Stats(ProjectName = OriginalSampleName, List_of_stats = Statistics_File_List)
    Final_Blast_File = Combine_Blastoutputs(ProjectName = OriginalSampleName, List_of_Blast_Files =BlastN_Output___List , Blast_Type = "blastN")
    Final_Fasta_File = Combine_Consensus_Fastas(ProjectName = OriginalSampleName, List_of_fastas = Consensus_Fasta_List)

    Move_File(Date = OriginalSampleName + "_RawFiles_", ProjectName = ProjectName, List_of_Files = Consensus_Fasta_List)
    Move_File(Date = OriginalSampleName + "_RawFiles_", ProjectName = ProjectName, List_of_Files = BlastN_Output___List)
    Move_File(Date = OriginalSampleName + "_RawFiles_", ProjectName = ProjectName, List_of_Files = Statistics_File_List)
    Move_File(Date = OriginalSampleName + "_RawFiles_", ProjectName = ProjectName, List_of_Files = Qual_Graph_File_List)
    Move_File(Date = OriginalSampleName + "_RawFiles_", ProjectName = ProjectName, List_of_Files = Cove_Graph_File_List)
    Move_File(Date = OriginalSampleName + "_RawFiles_", ProjectName = ProjectName, List_of_Files = BAM________File_List)
    Move_File(Date = OriginalSampleName + "_RawFiles_", ProjectName = ProjectName, List_of_Files = SNP_csv____File_List)

    Final_File_List = []
    Final_File_List = [Final_Graph_File_1, Final_Graph_File_2, Final_Blast_File, Final_Fasta_File]
    Move_File(Date = OriginalSampleName + "_SummaryFiles_", ProjectName = ProjectName, List_of_Files = Final_File_List)

    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Complete " + Directory + "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", flush = True)
    del(NewMetaDF)
    Sample_Counter += 1

print("Pipeline Complete", flush = True)