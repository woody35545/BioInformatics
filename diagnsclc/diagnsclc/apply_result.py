from tensorflow.python.keras.utils import np_utils
from tensorflow.python.keras.models import Sequential
from tensorflow.python.keras.layers import Dense, Activation
from tensorflow.python.keras.models import load_model
from itertools import product
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from itertools import product
from os.path import join
import tensorflow as tf
from tensorflow.keras import layers, models, optimizers
from tensorflow.python.keras.models import load_model




f = open('patient.txt', 'r')
test_patient = f.readline()

test_patient = K_mer(test_patient)

probability_model = tf.keras.Sequential([model, tf.keras.layers.Softmax()])
xhat = test_patient.reshape(-1,5597,125,1)

yhat = probability_model.predict(xhat)
print(yhat)


f = open('normal.txt', 'r')
test_normal = f.readline()

test_normal = K_mer(test_normal)

probability_model = tf.keras.Sequential([model, tf.keras.layers.Softmax()])
xhat = test_normal.reshape(-1,5597,125,1)

yhat = probability_model.predict(xhat)
print(yhat)


def apply(data):
   model = load_model('./save')
   encoded_data = encoder.K_mer(data)
   # encoded_data = encoder.onehot(data)
   probability_model = tf.keras.Sequential([model, tf.keras.layers.Softmax()])
   xhat = encoded_data.reshape(-1, 5597, 125, 1)

   result = probability_model.predict(xhat)
   print(result)

   if result[1]>0.7:
      print("You are patient")
      f = open('biomarker_info.txt', 'w')
      f.write(result[1])
      f.close()
   else:
      print("You are not patient")