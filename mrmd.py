#!/usr/bin/env python3
# -*- coding=utf-8 -*-

from sklearn.preprocessing import MinMaxScaler
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import LinearSVC
from feature_Rank import feature_rank
import argparse
import sklearn.metrics
import time
import logging
import os,sys
from math import ceil
from sklearn.preprocessing import LabelBinarizer

class Dim_Rd(object):
    def __init__(self, file_csv, args):
        self.file_csv = file_csv

        self.logger = args
    def read_data(self):  # default csv

        def read_csv():
            self.df = pd.read_csv(self.file_csv, engine='python').dropna(axis=1)
            datas = np.array(self.df)
            self.datas = datas
            self.X = datas[:, 1:]
            self.y = datas[:, 0]

        file_type = self.file_csv.split('.')[-1]
        if file_type == 'csv':
            read_csv()

    def range_steplen(self, start=1, end=1, length=1):
        self.start = start
        self.end = end
        self.length = length

    def Randomforest(self, X, y):

        clf = RandomForestClassifier(random_state=1, n_estimators=100,n_jobs=-1)
        ypred = sklearn.model_selection.cross_val_predict(clf, X, y, n_jobs=-1, cv=2)
        f1 = sklearn.metrics.f1_score(y, ypred, average='weighted')
        precision = sklearn.metrics.precision_score(self.y, ypred, average='weighted')
        recall = sklearn.metrics.recall_score(self.y, ypred, average='weighted')
        acc = sklearn.metrics.accuracy_score(self.y, ypred)
        lb = LabelBinarizer()
        lb.fit(self.y)

        y = lb.transform(self.y)
        ypred = lb.transform(ypred)
        auc = sklearn.metrics.roc_auc_score(y, ypred)

        return acc, f1, precision, recall, auc, ypred

    def Result(self, seqmax, clf, features, csvfile,args):
        ypred = sklearn.model_selection.cross_val_predict(clf, self.X[:, seqmax], self.y, n_jobs=-1, cv=5)

        cm = pd.crosstab(pd.Series(self.y, name='Actual'), pd.Series(ypred, name='Predicted'))
        f1 = sklearn.metrics.f1_score(self.y, ypred, average='weighted')
        acc = sklearn.metrics.accuracy_score(self.y, ypred, )
        precision = sklearn.metrics.precision_score(self.y, ypred, average='weighted')
        recall = sklearn.metrics.recall_score(self.y, ypred, average='weighted')
        lb = LabelBinarizer()
        y = lb.fit_transform(self.y)
        ypred = lb.transform(ypred)
        auc = sklearn.metrics.roc_auc_score(y, ypred)

        columns_index = [0]
        columns_index.extend([i + 1 for i in seqmax])
        data = np.concatenate((self.y.reshape(len(self.y), 1), self.X[:, seqmax]), axis=1)
        features_list = (self.df.columns.values)

        if args.n == -1:
            pass
        else:
            columns_index = columns_index[0:args.n + 1]
            data = data[:, 0:args.n + 1]
        df = pd.DataFrame(data, columns=features_list[columns_index])
        df.iloc[0, :].astype(int)
        df.to_csv(csvfile, index=None)

    def Dim_reduction(self, features, features_sorted, outfile, csvfile,args):
        features_number = []
        for value in features_sorted:
            f0=value[0]
            f=features[value[0]]
            f1 = features[value[0]]-1
            features_number.append(features[value[0]] - 1)
        stepSum = 0
        max = 0
        seqmax = []
        predmax = []
        scorecsv = outfile

        with open(scorecsv, 'w') as f:
            f.write('length,accuracy,f1,precision,recall,roc\n')
            for i in range(int(ceil((self.end - self.start) / self.length)) + 1):
                if (stepSum + self.start) < self.end:
                    n = stepSum + self.start
                else:
                    n = self.end

                stepSum += self.length

                ix = features_number[self.start - 1:n]
                acc, f1, precision, recall, auc, ypred = self.Randomforest(self.X[:, ix], self.y)
                benchmark = f1
                if benchmark > max:
                    max = benchmark
                    seqmax = ix

                f.write('{},{:0.4f},{:0.4f},{:0.4f},{:0.4f},{:0.4f}\n'.format(len(ix), acc, f1, precision, recall, auc))

        clf = RandomForestClassifier(random_state=1, n_estimators=100,n_jobs=-1)
        self.Result(seqmax, clf, features, csvfile,args)

    def run(self, inputfile,args):
        file = inputfile
        if args.m == None:
            args.m = ''.join(os.path.basename(args.i).split('.')[:-1]) + '.metrics.csv'

        metrics_file = args.m
        csvfile = args.o
        mrmr_featurLen = args.n
        features, features_sorted = feature_rank(file,  mrmr_featurLen,args.r)
        self.read_data()
        if int(args.e) == -1:
            args.e = len(pd.read_csv(file, engine='python').columns) - 1
        self.range_steplen(1, args.e, args.l)
        metrics_file = os.getcwd() + os.sep + 'Results' + os.sep + metrics_file
        csvfile = os.getcwd() + os.sep + 'Results' + os.sep + csvfile
        self.Dim_reduction(features, features_sorted, metrics_file, csvfile,args)


def mrmd(file,args):

    d = Dim_Rd(file, args)
    d.run(file,args)

