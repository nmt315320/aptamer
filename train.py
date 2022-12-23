from model.utils import *
from model.DL_ClassifierModel import *
import argparse
import numpy
from RicEnnTool import *
from iLearn.pubscripts import read_fasta_sequences, save_file
def DNA_to_array(DNA, module):
    DNA_array = []
    sum = 0
    for i in DNA:
        for j in i:
            sum = sum + module.wv[j]
        DNA_array.append(sum)
    return numpy.array(DNA_array)
def k_mergenerator(WaiTingKmer, k ,name):
    result = open("{}.txt".format(name), mode="w")
    DNA = []
    for i in WaiTingKmer:
        line = []
        for j in range(len(i)-k+1):
            line.append(i[j:j+k])
            result.write(i[j:j+k])
            result.write(" ")
        DNA.append(line)
        result.write("\n")
    return DNA


def read_fasta(fasta_file):
    try:
        fp = open(fasta_file)
    except IOError:
        print('cannot open ' + fasta_file + ', check if it exist!')
        exit()
    else:
        fp = open(fasta_file)
        lines = fp.readlines()

        fasta_dict = {}  # record seq for one id
        idlist = []  # record id list sorted
        gene_id = ""
        seq = ""
        s = []
        for line in lines:
            line = line.replace('\r', '')
            if line[0] == '>':
                if gene_id != "":
                    fasta_dict[gene_id] = seq.upper()
                    idlist.append(gene_id)

                gene_id = line.strip('\n')  # line.split('|')[1] all in > need to be id
            else:
                seq += line.strip('\n')
                s.append(seq)
        fasta_dict[gene_id] = seq.upper()  # last seq need to be record
        idlist.append(gene_id)

    return s, idlist
parser = argparse.ArgumentParser()
parser.add_argument('--k', default='6')
parser.add_argument('--d', default='32')
parser.add_argument('--s', default='64')
parser.add_argument('--f', default='128')
parser.add_argument('--metrics', default='MaF')
parser.add_argument('--device', default='cuda:0')
parser.add_argument('--savePath', default='out/')
args = parser.parse_args()


if __name__ == '__main__':
    k,d,s,f = int(args.k),int(args.d),int(args.s),int(args.f)
    device,path = args.device,args.savePath
    metrics = args.metrics
    full_path = os.path.realpath(__file__)
    curdir = os.path.dirname(full_path)
    modelPath = os.path.join(curdir, 'model')
    resultPath = os.path.join(curdir, 'result')
    dataPath = os.path.join(curdir, 'data')

    report = ["ACC", "MaF", "MiAUC", "MaAUC"]
    tezhenArray = ['DAC', 'DACC', 'DCC', 'TAC', 'TACC', 'TCC']
    fileTrainName = './data/dtrain.txt'
    fileTestName = './data/dtrain.txt'
    seq, id = read_fasta(fileTrainName)
    dataClass = DataClass('./data/data.txt', 0.2, 0.0, kmers=k)
    dataClass.vectorize("char2vec", feaSize=d, loadCache=True)

    # train_x1, train_y = xiangguanTezhenTiqu(os.path.join(dataPath, fileTrainName), tezhenArray)
    # test_x1, test_y = xiangguanTezhenTiqu(os.path.join(dataPath, fileTestName), tezhenArray)
    model = TextClassifier_SPPCNN(classNum=2, embedding=dataClass.vector['embedding'], SPPSize=s, feaSize=d, filterNum=f,
                                   contextSizeList=[1,3,5], hiddenList=[], 
                                   embDropout=0.3, fcDropout=0.5, useFocalLoss=True, weight=None, 
                                   )
    model.cv_train(dataClass, trainSize=1, batchSize=16, stopRounds=-1, earlyStop=10,
                   epoch=100, lr=0.0001, kFold=5, savePath=f'{path}CNN_s{s}_f{f}_k{k}_d{d}', report=report)