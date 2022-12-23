# encoding:utf-8
import os
from datetime import datetime

import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from tensorflow import keras

from tensorflow.keras.models import *
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l1, l2
from tensorflow.keras import backend as K

from tensorflow.keras import initializers
from tensorflow.keras import layers,Sequential,optimizers,datasets,metrics
from tensorflow.python.keras import Input
from tensorflow.python.keras.callbacks import Callback

from tensorflow.python.keras.layers import Embedding, Conv1D, MaxPooling1D, BatchNormalization, Dropout, \
    Bidirectional, Dense, GRU, Flatten


class AttLayer(layers.Layer):
    def __init__(self, attention_dim=50):
        # self.init = initializers.get('normal')
        self.init = initializers.RandomNormal(seed=10)
        self.supports_masking = True
        self.attention_dim = attention_dim
        super(AttLayer, self).__init__()

    def build(self, input_shape):
        assert len(input_shape) == 3
        self.W = K.variable(self.init((input_shape[-1], self.attention_dim)))
        self.b = K.variable(self.init((self.attention_dim, )))
        self.u = K.variable(self.init((self.attention_dim, 1)))
       # self.trainable_weights = [self.W, self.b, self.u]
        super(AttLayer, self).build(input_shape)

    def compute_mask(self, inputs, mask=None):
        return mask

    def call(self, x, mask=None):
        # size of x :[batch_size, sel_len, attention_dim]
        # size of u :[batch_size, attention_dim]
        # uit = tanh(xW+b)
        uit = K.tanh(K.bias_add(K.dot(x, self.W), self.b))
        ait = K.dot(uit, self.u)
        ait = K.squeeze(ait, -1)

        ait = K.exp(ait)

        if mask is not None:
            # Cast the mask to floatX to avoid float64 upcasting in theano
            ait *= K.cast(mask, K.floatx())
        ait /= K.cast(K.sum(ait, axis=1, keepdims=True) + K.epsilon(), K.floatx())
        ait = K.expand_dims(ait)
        weighted_input = x * ait
        output = K.sum(weighted_input, axis=1)

        return output

    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[-1])

    def get_config(self):
        config = {"attention_dim": self.attention_dim}#
        base_config = super(AttLayer, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))


class RicENNNWModel2(keras.Model):#inherit keras.Model
    def __init__(self):  # init
        super(RicENNNWModel2, self).__init__()

        # self.fc1=Embedding(num_words, 20, trainable=True)
        self.fc2 = Conv1D(filters=10, kernel_size=5, padding="valid", activation='relu')
        self.fc3 = MaxPooling1D(pool_size=4, strides=4)
        self.fc4 = BatchNormalization()

        self.fc5 = Dropout(0.5)
        self.fc6 = Bidirectional(GRU(50, return_sequences=True))
        self.fc7 = AttLayer(50)
        self.fc8 = Flatten()
        self.fc9 = Dense(2, activation='sigmoid')

    def call(self, inputs, training=None):  # compute inputs.shape = [b,28*28]
        # emb_en = self.fc1(inputs)
        enhancer_conv_layer = self.fc2(inputs)
        enhancer_max_pool_layer = self.fc3(enhancer_conv_layer)
        bn = self.fc4(enhancer_max_pool_layer)
        dt = self.fc5(bn)
        l_gru = self.fc6(dt)
        l_gru = tf.nn.tanh(l_gru) + l_gru
        l_att = self.fc7(l_gru)
        inner = self.fc8(l_att)
        preds = self.fc9(inner)
        return preds
