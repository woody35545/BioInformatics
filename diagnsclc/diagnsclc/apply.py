from tensorflow.python.keras.utils import np_utils
from tensorflow.python.keras.models import Sequential
from tensorflow.python.keras.layers import Dense, Activation
from tensorflow.python.keras.models import load_model
import encoder


def apply(data):
    
    model = load_model('./save')
    encoded_data = encoder.K_mer(data)
    # encoded_data = encoder.onehot(data)

    result = model.predict(encoded_data)
    print(result)