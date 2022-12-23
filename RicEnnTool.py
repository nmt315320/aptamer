# encoding:utf-8
# !/usr/bin/python
# _*_ coding:utf-8 _*_

import csv

import os
import sys
import re
from argparse import Namespace

import xlrd
#from Bio import SeqIO
import numpy as np
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve, auc, precision_recall_curve
import matplotlib as mlp
os.path.join(os.path.dirname(__file__), '../../')
sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from iLearn.descnucleotide import check_parameters, ACC
from iLearn.pubscripts import read_fasta_sequences, save_file

mlp.use('Agg')
import matplotlib.pyplot as plt

#from DNA import DNC
#from kmer import Kmer
from tensorflow import keras
import tensorflow as tf

def read_fa(fa):
	"""读入FASTA格式的文件对象fa，将其中的序列用字典结构存储"""
	seq = {}
	for line in fa:
		if line.startswith('>'):
			seqid = line.replace('>', '').split()[0]
			seq[seqid] = ''
		else:
			seq[seqid] += line.strip()
	return seq

def gc_content(s):
	"""读入DNA序列s，计算s中的GC含量"""
	gcc = (s.count('C')+s.count('G'))/len(s)*100
	return gcc

def read_tsv(fn):
    data = []
    f = open(fn, 'r')
    # f.readline()

    for line in f:
        sl = line.split()
        data.append(sl)

    return data


# function for string float conversion
def str_to_float_tezhen(arr):
    t = []
    for i in range(0, len(arr)):
        if arr[i]!='':
            t.append(float(arr[i]))

    return t


# function for mapping read data into features and target values
def get_map(data,index=2,dict={
          'Yes':[1.0, 0.0],
         'No':[0.0, 1.0]
    }):
    print("the representation of the labels")
    print(dict)

    mapped_data_x = []
    mapped_data_y = []
    mapped_data_bianhao=[]

    for d in data:
        mapped_data_bianhao_t=''
        mapped_data_x_t=''
        for d1 in d[index:]:
            d2=d1.split(':')
            arr=float(d2[1])
            #print("d2[1]:%s,arr:%0.3f"%(d2[1],arr))
            mapped_data_x_t=mapped_data_x_t+'%0.3f '%float(d2[1])
            #mapped_data_x.append('%0.3f'%float(d2[1]))#str_to_float(d2[1]))
            mapped_data_bianhao_t=mapped_data_bianhao_t+d2[0]+' '

        mapped_data_x.append(str_to_float_tezhen(mapped_data_x_t.split(' ')))
        mapped_data_bianhao.append(mapped_data_bianhao_t)
        mapped_data_y.append(dict[d[0]])

    return dict, mapped_data_x, mapped_data_y

'''
利用iLearn-master的工具 提取相关性特征
'''
def xiangguanTezhenTiqu(fileSrc,tezhenArray):
    args=Namespace();
    #提取各个特征
    args.all_index = False
    args.file = fileSrc
    args.format = 'svm'
    args.index = None
    args.lag = 5
    args.method =None# 'TACC'
    args.out = None#'E:/pythonwork/pythonwork/Bioinformatics/Enhancer/GaoEnhangcer/data/06/TACC/fileTest06_1.txt'
    args.type = 'DNA'
    args.udi = None
    fileSrcName,fileSrcExt= os.path.splitext(fileSrc)#分离文件名和文件扩展名
    fileSrcPath = os.path.dirname(fileSrc);
    for tezhen in tezhenArray:
        args.method=tezhen
        args.out=os.path.join(fileSrcPath,fileSrcName+"_"+tezhen+'.txt')
        print("args.out==%s"%args.out)
        my_property_name, my_property_value, kmer = check_parameters.check_acc_arguments(args)
        fastas = read_fasta_sequences.read_nucleotide_sequences(args.file)
        encodings = []
        if args.method == 'DAC' or args.method == 'TAC':
            encodings = ACC.make_ac_vector(fastas, my_property_name, my_property_value, args.lag, kmer)
        if args.method == 'DCC' or args.method == 'TCC':
            encodings = ACC.make_cc_vector(fastas, my_property_name, my_property_value, args.lag, kmer)
        if args.method == 'DACC' or args.method == 'TACC':
            encodings = ACC.make_acc_vector(fastas, my_property_name, my_property_value, args.lag, kmer)
        out_file = args.out if args.out != None else 'encoding.txt'
        save_file.save_file(encodings, args.format, out_file)

    test_x1 = None
    for index in range(0, len(tezhenArray)):
        tezhen = tezhenArray[index]
        fileTestName =  os.path.join(fileSrcPath, fileSrcName+"_"+tezhen+'.txt')
        di = {'Yes': [1.0, 0.0], 'No': [0.0, 1.0]}

        test_data = read_tsv(fileTestName)
        diTestleaf, test_x, test_y = get_map(test_data, index=1, dict=di)

        if index > 0:
            test_x1 = np.hstack((test_x1, np.array(test_x)))
        else:
            test_x1 = np.array(test_x)
    return  test_x1, np.array(test_y)
def canshu(yt, yp):
    tn, fp, fn, tp = confusion_matrix(yt, yp).ravel()
    acc = (tp + tn) / (tp + tn + fp + fn)
    sp = tn / (tn + fp)
    sn = tp / (tp + fn)
    print("tp=%0.3f,fn=%0.3f" % (tp, fn))

    pr = tp / (tp + fp)  # precision

    rc = tp / (tp +  fn)  # Recall
    npv=tn/(tn+fn)#NPV
    return acc, sp, sn, pr,npv #rc,npv


def zhuanhuan(pre, yuzhi):
    for evn in range(0, len(pre)):
        if pre[evn] <= yuzhi:
            pre[evn] = 0
        else:
            pre[evn] = 1



'''绘制ROC'''


def hzRoc(tuPath, fprArray, tprArray, roc_aucArray, colorArray, modelNameArray):
    font1 = {'family': 'DejaVu Sans',
             'weight': 'bold',
             'size': 14,
             }
    plt.figure()
    lw = 2
    plt.figure(figsize=(8, 8))
    for index in range(0, len(roc_aucArray)):
        plt.plot(fprArray[index], tprArray[index], color=colorArray[index],
                 lw=lw, label='%s(AUROC = %0.2f)' % (
                modelNameArray[index], roc_aucArray[index]))  ###recall为横坐标，precision为纵坐标做曲线
    color1 = 'lightgrey'
    plt.plot([0, 1], [0, 1], color=color1, lw=0.5, linestyle='--', label='Luck')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.02])
    plt.xlabel('False Positive Rate', font1)
    plt.ylabel('True Positive Rate', font1)
    # plt.title('The PRC curves of All Model on 5-mer Feature',font1)
    plt.legend(loc="lower right", prop=font1)
    ###设置坐标轴的粗细

    ax = plt.gca();  # 获得坐标轴的句柄
    ax.tick_params(labelsize=14)
    ax.spines['bottom'].set_linewidth(lw);  ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(lw);  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(lw);  ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(lw);  ####设置上部坐标轴的粗细
    ax.spines['bottom'].set_color(color1);  ###设置底部坐标轴的粗细
    ax.spines['left'].set_color(color1);  ####设置左边坐标轴的粗细
    ax.spines['right'].set_color(color1);  ###设置右边坐标轴的粗细
    ax.spines['top'].set_color(color1);  ####设置上部坐标轴的粗细
    plt.savefig(tuPath)
    plt.show()


'''绘制PRC'''


def hzPrc(tuPath, recallArray, precisionArray, prc_aucArray, colorArray, modelNameArray):
    font1 = {'family': 'DejaVu Sans',
             'weight': 'bold',
             'size': 14,
             }
    plt.figure()
    lw = 2
    plt.figure(figsize=(8, 8))
    for index in range(0, len(prc_aucArray)):
        plt.plot(recallArray[index], precisionArray[index], color=colorArray[index],
                 lw=lw, label='%s(AUPRC = %0.2f)' % (
                modelNameArray[index], prc_aucArray[index]))  ###recall为横坐标，precision为纵坐标做曲线
    color1 = 'lightgrey'
    plt.plot([0, 1], [1, 0], color=color1, lw=0.5, linestyle='--', label='Luck')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.02])

    plt.xlabel('Recall', font1)
    plt.ylabel('Precision', font1)
    # plt.title('The PRC curves of All Model on 5-mer Feature',font1)
    plt.legend(loc="lower left", prop=font1)
    ###设置坐标轴的粗细
    ax = plt.gca();  # 获得坐标轴的句柄
    ax.tick_params(labelsize=14)
    ax.spines['bottom'].set_linewidth(lw);  ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(lw);  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(lw);  ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(lw);  ####设置上部坐标轴的粗细
    ax.spines['bottom'].set_color(color1);  ###设置底部坐标轴的粗细
    ax.spines['left'].set_color(color1);  ####设置左边坐标轴的粗细
    ax.spines['right'].set_color(color1);  ###设置右边坐标轴的粗细
    ax.spines['top'].set_color(color1);  ####设置上部坐标轴的粗细
    plt.savefig(tuPath)
    plt.show()


'''绘制ROC'''


def hzRoc(tuPath, fprArray, tprArray, roc_aucArray, colorArray, modelNameArray):
    font1 = {'family': 'DejaVu Sans',
             'weight': 'bold',
             'size': 14,
             }
    plt.figure()
    lw = 2
    plt.figure(figsize=(8, 8))
    for index in range(0, len(roc_aucArray)):
        plt.plot(fprArray[index], tprArray[index], color=colorArray[index],
                 lw=lw, label='%s(AUROC = %0.2f)' % (
                modelNameArray[index], roc_aucArray[index]))  ###recall为横坐标，precision为纵坐标做曲线
    color1 = 'lightgrey'
    plt.plot([0, 1], [0, 1], color=color1, lw=0.5, linestyle='--', label='Luck')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.02])
    plt.xlabel('False Positive Rate', font1)
    plt.ylabel('True Positive Rate', font1)
    # plt.title('The PRC curves of All Model on 5-mer Feature',font1)
    plt.legend(loc="lower right", prop=font1)
    ###设置坐标轴的粗细

    ax = plt.gca();  # 获得坐标轴的句柄
    ax.tick_params(labelsize=14)
    ax.spines['bottom'].set_linewidth(lw);  ###设置底部坐标轴的粗细
    ax.spines['left'].set_linewidth(lw);  ####设置左边坐标轴的粗细
    ax.spines['right'].set_linewidth(lw);  ###设置右边坐标轴的粗细
    ax.spines['top'].set_linewidth(lw);  ####设置上部坐标轴的粗细
    ax.spines['bottom'].set_color(color1);  ###设置底部坐标轴的粗细
    ax.spines['left'].set_color(color1);  ####设置左边坐标轴的粗细
    ax.spines['right'].set_color(color1);  ###设置右边坐标轴的粗细
    ax.spines['top'].set_color(color1);  ####设置上部坐标轴的粗细
    plt.savefig(tuPath)
    plt.show()


'''绘制Prc和Roc'''


def hzPrcRoc(tuPath, fpr, tpr, roc_auc, coloritem, modelName, recall, precision, prc_auc):
    font1 = {'family': 'DejaVu Sans',
             'weight': 'bold',
             'size': 11,
             }
    plt.figure()
    lw = 2
    plt.figure(figsize=(10, 5))
    ax1 = plt.subplot(1, 2, 2)

    lwB = 0.5
    ax1.plot(fpr, tpr, color=coloritem,
             lw=lw, label='%s(AUROC = %0.3f)' % (
            modelName, roc_auc))  ###recall为横坐标，precision为纵坐标做曲线
    color1 = 'grey'
    ax1.plot([0, 1], [0, 1], color=color1, lw=lwB, label='Luck')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.02])
    # ax1.xlabel('False Positive Rate', font1)
    # ax1.ylabel('True Positive Rate', font1)
    # plt.title('The PRC curves of All Model on 5-mer Feature',font1)
    ax1.legend(loc="lower right", prop=font1)
    ###设置坐标轴的粗细

    # ax1 = plt.gca();  # 获得坐标轴的句柄
    ax1.tick_params(labelsize=14)
    ax1.spines['bottom'].set_linewidth(lwB);  ###设置底部坐标轴的粗细
    ax1.spines['left'].set_linewidth(lwB);  ####设置左边坐标轴的粗细
    ax1.spines['right'].set_linewidth(lwB);  ###设置右边坐标轴的粗细
    ax1.spines['top'].set_linewidth(lwB);  ####设置上部坐标轴的粗细
    ax1.spines['bottom'].set_color(color1);  ###设置底部坐标轴的粗细
    ax1.spines['left'].set_color(color1);  ####设置左边坐标轴的粗细
    ax1.spines['right'].set_color(color1);  ###设置右边坐标轴的粗细
    ax1.spines['top'].set_color(color1);  ####设置上部坐标轴的粗细

    ax2 = plt.subplot(1, 2, 1)
    ax2.plot(recall, precision, color=coloritem,
             lw=lw, label='%s(AUPRC = %0.3f)' % (
            modelName, prc_auc))  ###recall为横坐标，precision为纵坐标做曲线
    color1 = 'grey'
    ax2.plot([0, 1], [1, 0], color=color1, lw=lwB, label='Luck')

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.02])
    # ax2.xlabel('Recall', font1)
    # ax2.ylabel('Precision', font1)
    # plt.title('The PRC curves of All Model on 5-mer Feature',font1)
    ax2.legend(loc="lower left", prop=font1)
    ###设置坐标轴的粗细
    # ax2 = plt.gca();  # 获得坐标轴的句柄
    ax2.tick_params(labelsize=14)
    ax2.spines['bottom'].set_linewidth(lwB);  ###设置底部坐标轴的粗细
    ax2.spines['left'].set_linewidth(lwB);  ####设置左边坐标轴的粗细
    ax2.spines['right'].set_linewidth(lwB);  ###设置右边坐标轴的粗细
    ax2.spines['top'].set_linewidth(lwB);  ####设置上部坐标轴的粗细
    ax2.spines['bottom'].set_color(color1);  ###设置底部坐标轴的粗细
    ax2.spines['left'].set_color(color1);  ####设置左边坐标轴的粗细
    ax2.spines['right'].set_color(color1);  ###设置右边坐标轴的粗细
    ax2.spines['top'].set_color(color1);  ####设置上部坐标轴的粗细
    plt.savefig(tuPath)
    plt.show()


'''绘制Prc和Roc'''


def hzPrcRoc1(tuPath, fprArray, tprArray, roc_aucArray, colorArray, modelNameArray, recallArray, precisionArray,
              prc_aucArray):
    font1 = {'family': 'DejaVu Sans',
             'weight': 'bold',
             'size': 10,
             }
    plt.figure()
    lw = 2
    plt.figure(figsize=(12, 6))
    ax1 = plt.subplot(1, 2, 2)

    lwB = 0.5
    for index in range(0, len(roc_aucArray)):
        ax1.plot(fprArray[index], tprArray[index], color=colorArray[index],
                 lw=lw, label='%s(AUROC = %0.3f)' % (
                modelNameArray[index], roc_aucArray[index]))

    color1 = 'grey'
    ax1.plot([0, 1], [0, 1], color=color1, lw=lwB, label='Luck')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.02])
    # ax1.xlabel('False Positive Rate', font1)
    # ax1.ylabel('True Positive Rate', font1)
    # plt.title('The PRC curves of All Model on 5-mer Feature',font1)
    ax1.legend(loc="lower right", prop=font1)
    ###设置坐标轴的粗细

    # ax1 = plt.gca();  # 获得坐标轴的句柄
    ax1.tick_params(labelsize=14)
    ax1.spines['bottom'].set_linewidth(lwB);  ###设置底部坐标轴的粗细
    ax1.spines['left'].set_linewidth(lwB);  ####设置左边坐标轴的粗细
    ax1.spines['right'].set_linewidth(lwB);  ###设置右边坐标轴的粗细
    ax1.spines['top'].set_linewidth(lwB);  ####设置上部坐标轴的粗细
    ax1.spines['bottom'].set_color(color1);  ###设置底部坐标轴的粗细
    ax1.spines['left'].set_color(color1);  ####设置左边坐标轴的粗细
    ax1.spines['right'].set_color(color1);  ###设置右边坐标轴的粗细
    ax1.spines['top'].set_color(color1);  ####设置上部坐标轴的粗细

    ax2 = plt.subplot(1, 2, 1)
    for index in range(0, len(prc_aucArray)):
        ax2.plot(recallArray[index], precisionArray[index], color=colorArray[index],
                 lw=lw, label='%s(AUPRC = %0.3f)' % (
                modelNameArray[index], prc_aucArray[index]))

    color1 = 'grey'
    ax2.plot([0, 1], [1, 0], color=color1, lw=lwB, label='Luck')

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.02])
    # ax2.xlabel('Recall', font1)
    # ax2.ylabel('Precision', font1)
    # plt.title('The PRC curves of All Model on 5-mer Feature',font1)
    ax2.legend(loc="lower left", prop=font1)
    ###设置坐标轴的粗细
    # ax2 = plt.gca();  # 获得坐标轴的句柄
    ax2.tick_params(labelsize=14)
    ax2.spines['bottom'].set_linewidth(lwB);  ###设置底部坐标轴的粗细
    ax2.spines['left'].set_linewidth(lwB);  ####设置左边坐标轴的粗细
    ax2.spines['right'].set_linewidth(lwB);  ###设置右边坐标轴的粗细
    ax2.spines['top'].set_linewidth(lwB);  ####设置上部坐标轴的粗细
    ax2.spines['bottom'].set_color(color1);  ###设置底部坐标轴的粗细
    ax2.spines['left'].set_color(color1);  ####设置左边坐标轴的粗细
    ax2.spines['right'].set_color(color1);  ###设置右边坐标轴的粗细
    ax2.spines['top'].set_color(color1);  ####设置上部坐标轴的粗细
    plt.savefig(tuPath)
    plt.show()


'''
绘制柱状图
'''


def hzZzt(tuPath, yDataArray, xLabelArray, colorArray, modelNameArray):
    font1 = {'family': 'DejaVu Sans',
             'weight': 'bold',
             'size': 10,
             }

    lw = 2
    plt.figure(figsize=(20, 10))

    fig, ax = plt.subplots()

    x = np.arange(len(xLabelArray))
    plt.grid(color='lightgrey', axis='y')
    width = 0.2
    for index in range(0, len(
            yDataArray)):  # [(0.5707426376440461, 0.45902688860435337, 0.6824583866837388, 0.5578231292517006, 0.6824583866837388), (0.7693633414436334, 0.6518653690186537, 0.8868613138686131, 0.7181086849450008, 0.8868613138686131)]
        ax.bar(x + (index + 1) * width, yDataArray[index], width, alpha=0.9, label=modelNameArray[index],
               color=colorArray[index], zorder=2, hatch="---", ec='white')

    ax.set_xticks(x + (width * len(xLabelArray)) / 2)
    ax.set_xticklabels(xLabelArray, font1)
    # yLabelArray=['0.7','','0.8','','0.9','','1.0']
    # ax.set_yticklabels(yLabelArray, font1)

    ax.tick_params(length=0, labelsize=10)
    # plt.xlim([0.0, len(xLabelArray)])
    # plt.ylim([0.7, 1.02])
    plt.ylim([0.0, 1.0])
    plt.ylabel('Performance', font1)
    # plt.xlabel('Recall', font1)
    # plt.ylabel('Precision', font1)
    # plt.title('Enhancer prediction',font1)
    # plt.legend(loc="lower center", prop=font1)
    # 进行备注
    '''
    plt.annotate("Performance", (55, 20), xycoords='data',
                 xytext=(5, 38),
                 arrowprops=dict(arrowstyle='->'))
    '''
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)

    ax.spines['bottom'].set_color('grey')
    ax.spines['left'].set_visible(False)
    plt.legend(bbox_to_anchor=(0.0, -0.15, 1.0, -.15), loc=8, columnspacing=0,
               ncol=len(yDataArray), mode="expand", prop=font1, frameon=False, borderpad=0.5)
    plt.savefig(tuPath)
    plt.show()


'''绘制折线图'''


def hzZxt(tuPath, yDataArray, xLabelArray, colorArray, modelNameArray, tuNames=None):
    font1 = {'family': 'DejaVu Sans',
             'weight': 'bold',
             'size': 10,
             }
    lw = 2
    width = 0.2
    plt.figure(figsize=(10, 10))

    if tuNames != None and len(tuNames) == len(yDataArray):
        tuNum = len(tuNames)
        for numIndex in range(1, tuNum + 1):
            ax1 = plt.subplot(int(tuNum / 2 + 0.5), 2, int(numIndex % 2 + 1))
            # 有多少个model就有多少ACC
            x = np.arange(len(modelNameArray))
            for xianIndex in range(0, len(xLabelArray)):
                ax1.plot(x, yDataArray[numIndex - 1][xianIndex], label=xLabelArray[xianIndex])
                ax1.scatter(x, yDataArray[numIndex - 1][xianIndex], s=200)
            ax1.legend()
            # ax1.title(tuNames[numIndex-1], font1)
            ax1.set_xticks(x + width * (len(modelNameArray) + 1) / 2)
            ax1.set_xticklabels(modelNameArray, font1)

    # 绘制单个折线图
    elif (len(yDataArray) == 1):
        fig, ax = plt.subplots()
        x = np.arange(len(xLabelArray[0]))
        ax.grid(color='lightgrey', axis='y')
        # for xianIndex in range(0, len(xLabelArray)):
        ax.plot(x, yDataArray[0], label=xLabelArray[0])
        ax.scatter(x, yDataArray[0], s=200)
        # ax.set_xticks(x + width * (len(modelNameArray) + 1) / 2)
        ax.set_xticklabels(xLabelArray[0], font1)
        plt.title(modelNameArray)
        plt.xticks(x, xLabelArray[0])
        # ax.legend()

    # plt.ylim([0.0, 1.0])
    plt.ylabel('Performance', font1)
    # plt.xlabel('Recall', font1)
    # plt.ylabel('Precision', font1)
    # plt.title('Enhancer prediction',font1)
    # plt.legend(loc="lower center", prop=font1)
    # 进行备注
    '''
    plt.annotate("Performance", (55, 20), xycoords='data',
                 xytext=(5, 38),
                 arrowprops=dict(arrowstyle='->'))
    '''

    # plt.legend(bbox_to_anchor=(0.1, -0.15, 0.8, -.15), loc=3,columnspacing=0,
    #          ncol=len(yDataArray), mode="expand",  prop=font1,frameon=False,borderpad=0.5)
    plt.savefig(tuPath)
    plt.show()


# 获取工具预测结果的真实值和预测值
'''
iEnhancer_CNN分析结果统计
fileExcel:记录工具分析结果的excel文件
ytIndex:真实值所在列
ypIndex:预测值所在列
qshIndex:数据起始行，默认为1，即第二行开始
ytNe:真实值正样本标志
ypNe:预测值负样本标志
dict:字典
'''
def hqgjTp(fileExcel, ytIndex=0, ypIndex=1, qshIndex=1, ytNe='ne', ypNe='not', dict={'Yes': [1.0], 'No': [0.0]}):
    # 打开文件
    data = xlrd.open_workbook(fileExcel)
    table = data.sheet_by_index(0)
    yt = []  # 真实值
    yp = []  # 预测值
    for hang in range(qshIndex, table.nrows):
        ytStr = str(table.cell(hang, ytIndex).value)
        # 如果找到了负样本标志
        if len(re.findall(ytNe, ytStr, flags=re.IGNORECASE)) > 0:
            yt.append(dict['No'])
        else:
            yt.append(dict['Yes'])
        ypStr = str(table.cell(hang, ypIndex).value)
        if len(re.findall(ypNe, ypStr, flags=re.IGNORECASE)) > 0:
            yp.append(dict['No'])
        else:
            yp.append(dict['Yes'])
    return yt, yp


def iEnhancer_ELTp(dirPath, dict={'Yes': [1.0], 'No': [0.0]}, flag=None):
    yt = []  # 真实值
    yp = []  # 预测值
    if (flag == None):
        flag = '--------------------------------------------------------------------------------------------------------------------------------------------------------------------'
    for root, dirs, files in os.walk(dirPath):
        for file in files:
            # 获取文件所属目录
            # print(root)
            # 获取文件路径
            # print(os.path.join(root, file))
            fileIn = open(os.path.join(root, file), 'r')
            neirong = fileIn.read()
            neirong = neirong.split(flag)
            for index in range(0, len(neirong)):  # acc---cba
                neirongtemp = neirong[index]
                yangben = None
                yangbenJie = None
                if '_po_' in neirongtemp:
                    yangben = 'Yes'
                    yt.append(dict['Yes'])
                elif '_ne_' in neirongtemp:
                    yangben = 'No'
                    yt.append(dict['No'])
                if 'non enhancer' in neirongtemp:
                    yangbenJie = 'No'
                    yp.append(dict['No'])
                elif 'weak enhancer' in neirongtemp:
                    yangbenJie = 'Yes'
                    yp.append(dict['Yes'])
                elif 'strong enhancer' in neirongtemp:
                    yp.append(dict['Yes'])
                    yangbenJie = 'Yes'
                # print(neirongtemp)
                # print('%s判为%s'%(yangben,yangbenJie))
    return yt, yp


'''计算模型在不同的数据表示下的参数值'''


def duibiShuju(fileTrainName, fileTestName, fileValName):
    numberwords = 5187
    di = {'Yes': [1.0], 'No': [0.0]}
    train_data = read_tsv(fileTrainName)
    diTrain, train_x, train_y = get_map(train_data, index=1, dict=di)

    train_y = np.array(train_y)
    train_y_f = train_y

    test_data = read_tsv(fileTestName)
    diTestleaf, test_x, test_y = get_map(test_data, index=1, dict=di)
    test_y = np.array(test_y)

    test_y_f = test_y

    val_data = read_tsv(fileValName)
    diTestFlower, val_x, val_y = get_map(val_data, index=1, dict=di)
    val_y = np.array(val_y)

    val_y_f = val_y

    maxlen = 2188
    x_train = keras.preprocessing.sequence.pad_sequences(train_x, maxlen, padding='post')
    x_val = keras.preprocessing.sequence.pad_sequences(val_x, maxlen, padding='post')
    x_test = keras.preprocessing.sequence.pad_sequences(test_x, maxlen, padding='post')

    roc_auc_train_Array = []
    roc_auc_test_Array = []
    roc_auc_val_Array = []
    prc_auc_train_Array = []
    prc_auc_test_Array = []
    prc_auc_val_Array = []
    fpr_train_Array = []
    tpr_train_Array = []
    fpr_test_Array = []
    tpr_test_Array = []
    fpr_val_Array = []
    tpr_val_Array = []

    precision_train_Array = []
    recall_train_Array = []
    precision_test_Array = []
    recall_test_Array = []
    precision_val_Array = []
    recall_val_Array = []
    cansh_train_Array = []
    cansh_test_Array = []
    cansh_val_Array = []
    modelNames = []
    '''OneHotMer3Model3'''
    print('OneHotMer3Model3')
    modelNames.append('OM33')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer3Model3'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_train = predit_train[:, 0]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)
    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    '''OneHotMer4Model3'''

    print('OneHotMer4Model3')
    modelNames.append('OM43')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer4Model3'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]

    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)
    '''OneHotMer5Model3'''
    print('OneHotMer5Model3')
    modelNames.append('OM53')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer5Model3'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    '''OneHotMer345Model3'''
    print('OneHotMer345Model3')
    modelNames.append('OM3453')

    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer345Model3'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    '''OneHotMer3Model2'''
    print('OneHotMer3Model2')
    modelNames.append('OM32')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer3Model2'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_train = predit_train[:, 0]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)
    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    '''OneHotMer4Model2'''

    print('OneHotMer4Model2')
    modelNames.append('OM42')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer4Model2'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]

    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    '''OneHotMer5Model2'''
    print('OneHotMer5Model2')
    modelNames.append('OM52')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer5Model2'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    # 信息

    '''OneHotMer345Model2'''
    print('OneHotMer345Model2')
    modelNames.append('OM3452')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer345Model2'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    maxlen = 6555
    x_train = keras.preprocessing.sequence.pad_sequences(train_x, maxlen, padding='post')
    x_val = keras.preprocessing.sequence.pad_sequences(val_x, maxlen, padding='post')
    x_test = keras.preprocessing.sequence.pad_sequences(test_x, maxlen, padding='post')
    '''OneHotMer345hModel3'''
    print('OneHotMer345hModel3')
    modelNames.append('OM345h3')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer345hModel3'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    '''OneHotMer345hModel2'''
    print('OneHotMer345hModel2')
    modelNames.append('OM345h2')

    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer345hModel2'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    return roc_auc_train_Array, \
           roc_auc_test_Array, \
           roc_auc_val_Array, \
           prc_auc_train_Array, \
           prc_auc_test_Array, \
           prc_auc_val_Array, \
           fpr_train_Array, \
           tpr_train_Array, fpr_test_Array, \
           tpr_test_Array, \
           fpr_val_Array, \
           tpr_val_Array, \
           precision_train_Array, \
           recall_train_Array, \
           precision_test_Array, \
           recall_test_Array, \
           precision_val_Array, \
           recall_val_Array, \
           cansh_train_Array, \
           cansh_test_Array, \
           cansh_val_Array, \
           modelNames


'''计算单个模型的指标参数值,在验证集精度最好的观察点下的表现'''


def duibiShujuDange_valacc(fileTrainName, fileTestName, fileValName, modelName, checkpointpath, model):
    di = {'Yes': [1.0, 0.0], 'No': [0.0, 1.0]}
    train_data = read_tsv(fileTrainName)
    diTrain, train_x, train_y = get_map(train_data, index=1, dict=di)

    train_y = np.array(train_y)
    train_y_f = train_y[:, 0]

    test_data = read_tsv(fileTestName)
    diTestleaf, test_x, test_y = get_map(test_data, index=1, dict=di)
    test_y = np.array(test_y)

    test_y_f = test_y[:, 0]

    val_data = read_tsv(fileValName)
    diTestFlower, val_x, val_y = get_map(val_data, index=1, dict=di)
    val_y = np.array(val_y)

    val_y_f = val_y[:, 0]

    maxlen = 2188
    x_train = keras.preprocessing.sequence.pad_sequences(train_x, maxlen, padding='post')
    x_val = keras.preprocessing.sequence.pad_sequences(val_x, maxlen, padding='post')
    x_test = keras.preprocessing.sequence.pad_sequences(test_x, maxlen, padding='post')

    roc_auc_train_Array = []
    roc_auc_test_Array = []
    roc_auc_val_Array = []
    prc_auc_train_Array = []
    prc_auc_test_Array = []
    prc_auc_val_Array = []
    fpr_train_Array = []
    tpr_train_Array = []
    fpr_test_Array = []
    tpr_test_Array = []
    fpr_val_Array = []
    tpr_val_Array = []

    precision_train_Array = []
    recall_train_Array = []
    precision_test_Array = []
    recall_test_Array = []
    precision_val_Array = []
    recall_val_Array = []
    cansh_train_Array = []
    cansh_test_Array = []
    cansh_val_Array = []

    print(checkpointpath)

    # checkpoint = tf.train.Checkpoint(model=model)

    # checkpoint.restore(tf.train.load_checkpoint(checkpointpath))

    model.load_weights(checkpointpath)

    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)
    return roc_auc_train_Array, \
           roc_auc_test_Array, \
           roc_auc_val_Array, \
           prc_auc_train_Array, \
           prc_auc_test_Array, \
           prc_auc_val_Array, \
           fpr_train_Array, \
           tpr_train_Array, fpr_test_Array, \
           tpr_test_Array, \
           fpr_val_Array, \
           tpr_val_Array, \
           precision_train_Array, \
           recall_train_Array, \
           precision_test_Array, \
           recall_test_Array, \
           precision_val_Array, \
           recall_val_Array, \
           cansh_train_Array, \
           cansh_test_Array, \
           cansh_val_Array, \
           modelName


'''计算单个模型的指标参数值,在验证集精度最好的观察点下的表现'''


def duibiShujuDange(fileTrainName, fileTestName, fileValName, modelName, checkpointpath):
    di = {'Yes': [1.0, 0.0], 'No': [0.0, 1.0]}
    train_data = read_tsv(fileTrainName)
    diTrain, train_x, train_y = get_map(train_data, index=1, dict=di)

    train_y = np.array(train_y)
    train_y_f = train_y[:, 0]

    test_data = read_tsv(fileTestName)
    diTestleaf, test_x, test_y = get_map(test_data, index=1, dict=di)
    test_y = np.array(test_y)

    test_y_f = test_y[:, 0]

    val_data = read_tsv(fileValName)
    diTestFlower, val_x, val_y = get_map(val_data, index=1, dict=di)
    val_y = np.array(val_y)

    val_y_f = val_y[:, 0]

    maxlen = 2188
    x_train = keras.preprocessing.sequence.pad_sequences(train_x, maxlen, padding='post')
    x_val = keras.preprocessing.sequence.pad_sequences(val_x, maxlen, padding='post')
    x_test = keras.preprocessing.sequence.pad_sequences(test_x, maxlen, padding='post')

    roc_auc_train_Array = []
    roc_auc_test_Array = []
    roc_auc_val_Array = []
    prc_auc_train_Array = []
    prc_auc_test_Array = []
    prc_auc_val_Array = []
    fpr_train_Array = []
    tpr_train_Array = []
    fpr_test_Array = []
    tpr_test_Array = []
    fpr_val_Array = []
    tpr_val_Array = []

    precision_train_Array = []
    recall_train_Array = []
    precision_test_Array = []
    recall_test_Array = []
    precision_val_Array = []
    recall_val_Array = []
    cansh_train_Array = []
    cansh_test_Array = []
    cansh_val_Array = []

    print(modelName)

    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)
    return roc_auc_train_Array, \
           roc_auc_test_Array, \
           roc_auc_val_Array, \
           prc_auc_train_Array, \
           prc_auc_test_Array, \
           prc_auc_val_Array, \
           fpr_train_Array, \
           tpr_train_Array, fpr_test_Array, \
           tpr_test_Array, \
           fpr_val_Array, \
           tpr_val_Array, \
           precision_train_Array, \
           recall_train_Array, \
           precision_test_Array, \
           recall_test_Array, \
           precision_val_Array, \
           recall_val_Array, \
           cansh_train_Array, \
           cansh_test_Array, \
           cansh_val_Array, \
           modelName


'''计算模型在不同的数据表示下的参数值,这里只计算Mer345即345Mer竖向拼接后训练的CNN+RNN+ATT(Model2)
，去ATT(Model3),去CNN(Model4),去RNN（Model5）'''


def duibiShujuMer345(fileTrainName, fileTestName, fileValName):
    di = {'Yes': [1.0], 'No': [0.0]}
    train_data = read_tsv(fileTrainName)
    diTrain, train_x, train_y = get_map(train_data, index=1, dict=di)

    train_y = np.array(train_y)
    train_y_f = train_y

    test_data = read_tsv(fileTestName)
    diTestleaf, test_x, test_y = get_map(test_data, index=1, dict=di)
    test_y = np.array(test_y)

    test_y_f = test_y

    val_data = read_tsv(fileValName)
    diTestFlower, val_x, val_y = get_map(val_data, index=1, dict=di)
    val_y = np.array(val_y)

    val_y_f = val_y

    maxlen = 2188
    x_train = keras.preprocessing.sequence.pad_sequences(train_x, maxlen, padding='post')
    x_val = keras.preprocessing.sequence.pad_sequences(val_x, maxlen, padding='post')
    x_test = keras.preprocessing.sequence.pad_sequences(test_x, maxlen, padding='post')

    roc_auc_train_Array = []
    roc_auc_test_Array = []
    roc_auc_val_Array = []
    prc_auc_train_Array = []
    prc_auc_test_Array = []
    prc_auc_val_Array = []
    fpr_train_Array = []
    tpr_train_Array = []
    fpr_test_Array = []
    tpr_test_Array = []
    fpr_val_Array = []
    tpr_val_Array = []

    precision_train_Array = []
    recall_train_Array = []
    precision_test_Array = []
    recall_test_Array = []
    precision_val_Array = []
    recall_val_Array = []
    cansh_train_Array = []
    cansh_test_Array = []
    cansh_val_Array = []
    modelNames = []

    '''OneHotMer345Model2'''
    print('OneHotMer345Model2')
    modelNames.append('OM3452')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer345Model2'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    '''OneHotMer345Model3'''
    print('OneHotMer345Model3')
    modelNames.append('OM3453')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer345Model3'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    '''OneHotMer345Model4'''
    print('OneHotMer345Model4')
    modelNames.append('OM3454')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer345Model4'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)

    '''OneHotMer345Model5'''
    print('OneHotMer345Model5')
    modelNames.append('OM3455')
    checkpointpath = 'E:/pythonwork/Bioinformatics/Enhancer/myEnhancer/mymodel/./OneHotMer345Model5'
    model = keras.models.load_model(checkpointpath)
    predit_train = model.predict(x_train)
    predit_test = model.predict(x_test)
    predit_test = predit_test[:, 0]
    predit_train = predit_train[:, 0]
    test_y = np.array(test_y_f)  # test_y[:,1]
    train_y = np.array(train_y_f)  # train_y[:,1]

    predit_val = model.predict(x_val)
    predit_val = predit_val[:, 0]
    val_y = np.array(val_y_f)  # val_y[:,1]
    # 计算AUROC
    (fpr_train, tpr_train, thresholds_train) = roc_curve(train_y, predit_train, pos_label=[1.0])
    roc_auc_train = auc(fpr_train, tpr_train)
    roc_auc_train_Array.append(roc_auc_train)
    fpr_train_Array.append(fpr_train)
    tpr_train_Array.append(tpr_train)

    (fpr_test, tpr_test, thresholds_test) = roc_curve(test_y, predit_test, pos_label=[1.0])
    roc_auc_test = auc(fpr_test, tpr_test)
    roc_auc_test_Array.append(roc_auc_test)
    fpr_test_Array.append(fpr_test)
    tpr_test_Array.append(tpr_test)

    (fpr_val, tpr_val, thresholds_val) = roc_curve(val_y, predit_val, pos_label=[1.0])
    roc_auc_val = auc(fpr_val, tpr_val)
    roc_auc_val_Array.append(roc_auc_val)
    fpr_val_Array.append(fpr_val)
    tpr_val_Array.append(tpr_val)

    # 计算AUPRC
    precision_train, recall_train, thresholds = precision_recall_curve(train_y, predit_train, pos_label=[1.0])
    prc_auc_train = auc(recall_train, precision_train)  # 训练集集
    prc_auc_train_Array.append(prc_auc_train)
    precision_train_Array.append(precision_train)
    recall_train_Array.append(recall_train)

    precision_test, recall_test, thresholds = precision_recall_curve(test_y, predit_test, pos_label=[1.0])
    prc_auc_test = auc(recall_test, precision_test)  # 训练集集
    prc_auc_test_Array.append(prc_auc_test)
    precision_test_Array.append(precision_test)
    recall_test_Array.append(recall_test)

    precision_val, recall_val, thresholds = precision_recall_curve(val_y, predit_val, pos_label=[1.0])
    prc_auc_val = auc(recall_val, precision_val)  # 训练集集
    prc_auc_val_Array.append(prc_auc_val)
    precision_val_Array.append(precision_val)
    recall_val_Array.append(recall_val)

    zhuanhuan(predit_train, 0.5)
    zhuanhuan(predit_test, 0.5)
    zhuanhuan(predit_val, 0.5)

    canshuRe = canshu(test_y, predit_test)

    cansh_test_Array.append(canshuRe)
    canshuRe = canshu(train_y, predit_train)

    cansh_train_Array.append(canshuRe)

    canshuRe = canshu(val_y, predit_val)

    cansh_val_Array.append(canshuRe)
    return roc_auc_train_Array, \
           roc_auc_test_Array, \
           roc_auc_val_Array, \
           prc_auc_train_Array, \
           prc_auc_test_Array, \
           prc_auc_val_Array, \
           fpr_train_Array, \
           tpr_train_Array, fpr_test_Array, \
           tpr_test_Array, \
           fpr_val_Array, \
           tpr_val_Array, \
           precision_train_Array, \
           recall_train_Array, \
           precision_test_Array, \
           recall_test_Array, \
           precision_val_Array, \
           recall_val_Array, \
           cansh_train_Array, \
           cansh_test_Array, \
           cansh_val_Array, \
           modelNames



'''绘制训练过程的accuracy和loss曲线'''


def hzXunlianGc(history1, modelName, tuPath):
    print(history1.history.keys())
    # summarize history for accuracy
    ax1 = plt.subplot(1, 2, 1)
    plt.plot(history1.history['accuracy'])
    plt.plot(history1.history['val_accuracy'])
    plt.title(modelName + 'accuracy')
    plt.ylabel('accuracy')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')

    ax2 = plt.subplot(1, 2, 2)
    # summarize history for loss
    plt.plot(history1.history['loss'])
    plt.plot(history1.history['val_loss'])
    plt.title(modelName + 'loss')
    plt.ylabel('loss')
    plt.xlabel('epoch')
    plt.legend(['train', 'val'], loc='upper left')
    plt.savefig(tuPath)
    #plt.show()


