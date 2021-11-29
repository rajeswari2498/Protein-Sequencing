"""
Protein Sequencing Project
Name:
Roll Number:
"""

import hw6_protein_tests as test

project = "Protein" # don't edit this

### WEEK 1 ###

'''
readFile(filename)
#1 [Check6-1]
Parameters: str
Returns: str
'''
def readFile(filename):
    text=open(filename,"r")
    read_file=text.read()
    text.close()
    result=read_file.replace("\n","")
    return result


'''
dnaToRna(dna, startIndex)
#2 [Check6-1]
Parameters: str ; int
Returns: list of strs
'''
def dnaToRna(dna, startIndex):
    list=[]
    list1=[]
    for i in range(startIndex,len(dna),3):
        list.append(dna[i:i+3])
        if dna[i:i+3]=="TAG" or dna[i:i+3]=="TAA" or dna[i:i+3]=="TGA":
            break
    for each in list:
        each=each.replace("T","U")
        list1.append(each)
    # print(list1)
    return list1     


'''
makeCodonDictionary(filename)
#3 [Check6-1]
Parameters: str
Returns: dict mapping strs to strs
'''
def makeCodonDictionary(filename):  
    import json
    f=open(filename,"r")
    data=json.loads(f.read())
    # print(data)
    dictionary={}
    for key,value in data.items():
        for each in value:
            replacing=each.replace("T","U")
            dictionary[replacing]=key
    # print(dictionary)
    return dictionary


'''
generateProtein(codons, codonD)
#4 [Check6-1]
Parameters: list of strs ; dict mapping strs to strs
Returns: list of strs
'''
def generateProtein(codons, codonD):
    protein_list=[]
    for each in codons:
        if each=="AUG" and "Start" not in protein_list:
            protein_list.append("Start")
        else:
            protein_list.append(codonD[each])
    # print(protein_list) 
    return protein_list


'''
synthesizeProteins(dnaFilename, codonFilename)
#5 [Check6-1]
Parameters: str ; str
Returns: 2D list of strs
'''
def synthesizeProteins(dnaFilename, codonFilename):
    dna_read=readFile(dnaFilename)
    codon_dictionary=makeCodonDictionary(codonFilename)
    # print(dna_read,"\n",codon_dictionary)
    result_list=[]
    count=0
    j=0
    while j<len(dna_read):
        if dna_read[j:j+3]=="ATG":
            startIndex=j
            make_dnatorna=dnaToRna(dna_read,startIndex)
            make_rnatoprotein=generateProtein(make_dnatorna,codon_dictionary)
            result_list.append(make_rnatoprotein)
            j=j+(3*len(make_dnatorna))
        else:
            j=j+1
            count=count+1
    # print(count)
    # print(result_list)
    return result_list


def runWeek1():
    print("Human DNA")
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    print("Elephant DNA")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")


### WEEK 2 ###

'''
commonProteins(proteinList1, proteinList2)
#1 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs
Returns: 2D list of strs
'''
def commonProteins(proteinList1, proteinList2):
    list_protein=[]
    for each in proteinList1:
        if each in proteinList2 and each not in list_protein:
            list_protein.append(each)
    # print(list_protein) 
    return list_protein


'''
combineProteins(proteinList)
#2 [Check6-2]
Parameters: 2D list of strs
Returns: list of strs
'''
def combineProteins(proteinList):
    list_aminoacids=[]
    for each in proteinList:
        for i in each:
            list_aminoacids.append(i)
    # print(list_aminoacids)
    return list_aminoacids


'''
aminoAcidDictionary(aaList)
#3 [Check6-2]
Parameters: list of strs
Returns: dict mapping strs to ints
'''
def aminoAcidDictionary(aaList):
    dictionary={}
    for each in aaList:
        if each not in dictionary:
            dictionary[each]=aaList.count(each)
    # print(dictionary)
    return dictionary


'''
findAminoAcidDifferences(proteinList1, proteinList2, cutoff)
#4 [Check6-2]
Parameters: 2D list of strs ; 2D list of strs ; float
Returns: 2D list of values
'''
def findAminoAcidDifferences(proteinList1, proteinList2, cutoff):
    three_element_list=[]
    list1_protein=combineProteins(proteinList1)
    list1_amino_acid=aminoAcidDictionary(list1_protein)
    # print(proteinList1)
    list2_protein=combineProteins(proteinList2)
    list2_amino_acid=aminoAcidDictionary(list2_protein)
    # print(proteinList2)
    # print(cutoff)
    list1_diff=list(set(list2_protein) - set(list1_protein))
    list2_diff=list(set(list1_protein) - set(list2_protein))
    for each_diff1 in list1_diff:
        list1_amino_acid[each_diff1]=0
    for each_diff2 in list2_diff:
        list2_amino_acid[each_diff2]=0
    length1=len(list1_protein)
    aminoacid_freq1={}
    for each_1 in list1_amino_acid:
        aminoacid_freq1[each_1]=list1_amino_acid[each_1]/length1
    length2=len(list2_protein)
    aminoacid_freq2={}
    for each_2 in list2_amino_acid:
        aminoacid_freq2[each_2]=list2_amino_acid[each_2]/length2
    for key,value in aminoacid_freq1.items():
        list_aminoacid=[]
        if key!="Start" and key!="Stop":
            if key in aminoacid_freq2.keys():
                diff_frequencies=abs(aminoacid_freq1[key]-aminoacid_freq2[key])
                if diff_frequencies>cutoff:
                    list_aminoacid.append(key)
                    list_aminoacid.append(aminoacid_freq1[key])
                    list_aminoacid.append(aminoacid_freq2[key])
                    three_element_list.append(list_aminoacid) 
    # print(three_element_list)               
    return three_element_list


'''
displayTextResults(commonalities, differences)
#5 [Check6-2]
Parameters: 2D list of strs ; 2D list of values
Returns: None
'''
def displayTextResults(commonalities, differences):
    # print(commonalities)
    print("\n","These are the common proteins","\n")
    for each in sorted(commonalities):
        # print(each)
        temp=""
        for i in each:
            if i!="Start" and i!="Stop":
                temp=temp+"-"+i
        print(temp.strip("-"),"\n") 
    print("\n","These are the amino acids that occured at the most different rates","\n")
    # print(differences)  
    for percentage in differences:
        list1_percentage=percentage[1]
        list2_percentage=percentage[2]
        percentage_1="{:.2%}".format(list1_percentage)
        percentage_2="{:.2%}".format(list2_percentage)
        print(percentage[0],percentage_1,"in seq1,",percentage_2,"in seq2")
    return


def runWeek2():
    humanProteins = synthesizeProteins("data/human_p53.txt", "data/codon_table.json")
    elephantProteins = synthesizeProteins("data/elephant_p53.txt", "data/codon_table.json")

    commonalities = commonProteins(humanProteins, elephantProteins)
    differences = findAminoAcidDifferences(humanProteins, elephantProteins, 0.005)
    displayTextResults(commonalities, differences)


### WEEK 3 ###

'''
makeAminoAcidLabels(proteinList1, proteinList2)
#2 [Hw6]
Parameters: 2D list of strs ; 2D list of strs
Returns: list of strs
'''
def makeAminoAcidLabels(proteinList1, proteinList2):
    return


'''
setupChartData(labels, proteinList)
#3 [Hw6]
Parameters: list of strs ; 2D list of strs
Returns: list of floats
'''
def setupChartData(labels, proteinList):
    return


'''
createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None)
#4 [Hw6] & #5 [Hw6]
Parameters: list of strs ; list of floats ; str ; list of floats ; str ; [optional] list of strs
Returns: None
'''
def createChart(xLabels, freqList1, label1, freqList2, label2, edgeList=None):
    import matplotlib.pyplot as plt
    return


'''
makeEdgeList(labels, biggestDiffs)
#5 [Hw6]
Parameters: list of strs ; 2D list of values
Returns: list of strs
'''
def makeEdgeList(labels, biggestDiffs):
    return


'''
runFullProgram()
#6 [Hw6]
Parameters: no parameters
Returns: None
'''
def runFullProgram():
    return


### RUN CODE ###

# This code runs the test cases to check your work
if __name__ == "__main__":
    # test.testReadFile()
    # test.testDnaToRna()
    # test.testMakeCodonDictionary()
    # test.testGenerateProtein()
    # test.testSynthesizeProteins()
    # print("\n" + "#"*15 + " WEEK 1 TESTS " +  "#" * 16 + "\n")
    # test.week1Tests()
    # print("\n" + "#"*15 + " WEEK 1 OUTPUT " + "#" * 15 + "\n")
    # runWeek1()


    test.testCommonProteins()
    test.testCombineProteins()
    test.testAminoAcidDictionary()
    test.testFindAminoAcidDifferences()
    # Uncomment these for Week 2 ##
    print("\n" + "#"*15 + " WEEK 2 TESTS " +  "#" * 16 + "\n")
    test.week2Tests()
    print("\n" + "#"*15 + " WEEK 2 OUTPUT " + "#" * 15 + "\n")
    runWeek2()

    ## Uncomment these for Week 3 ##
    """
    print("\n" + "#"*15 + " WEEK 3 TESTS " +  "#" * 16 + "\n")
    test.week3Tests()
    print("\n" + "#"*15 + " WEEK 3 OUTPUT " + "#" * 15 + "\n")
    runFullProgram()
    """
