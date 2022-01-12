import pickle
import numpy as np
from scipy.io import loadmat
from tensorflow.keras.models import load_model
#import sys
from sklearn.preprocessing import StandardScaler
import scipy.io

DNN_Model = 30
#Ch = sys.argv[1]
#mod = sys.argv[2]
#scheme = sys.argv[3]
#hl1 = sys.argv[4]
#hl2 = sys.argv[5]
#hl3 = sys.argv[6]

Ch = 'VTV_UC' #sys.argv[1]
mod = 'QPSK' #sys.argv[2]
scheme = 'i11bd_noW' #sys.argv[3]
hl1 = 35 #sys.argv[6]
hl2 = 18 #sys.argv[7]
hl3 = 35 #sys.argv[8]

SNR_index = np.arange(1, 9)
j = 1

for j in SNR_index:
#if j == 1:
    mat = loadmat('./{}_DNN_Dataset_fd_1600_snr_{}.mat'.format(scheme,j))
    Testing_Dataset = mat['DNN_Datasets']
    Testing_Dataset = Testing_Dataset[0, 0]
    X = Testing_Dataset['Test_X']
    Y = Testing_Dataset['Test_Y']
    print('Loaded Dataset Inputs: ', X.shape)
    print('Loaded Dataset Outputs: ', Y.shape)
    # Normalizing Datasets
    scalerx = StandardScaler()
    scalerx.fit(X)
    scalery = StandardScaler()
    scalery.fit(Y)
    XS = scalerx.transform(X)
    YS = scalery.transform(Y)
    XS = XS.transpose()
    YS = YS.transpose()

    model = load_model('./{}_DNN_{}{}{}_{}.h5'.format(scheme, hl1, hl2, hl3, DNN_Model))
    print('Model Loaded: ', DNN_Model)
    # Testing the model

    Y_pred = model.predict(XS)

    XS = XS.transpose()
    YS = YS.transpose()
    Y_pred = Y_pred.transpose()

    Original_Testing_X = scalerx.inverse_transform(XS)
    Original_Testing_Y = scalery.inverse_transform(YS)
    Prediction_Y = scalery.inverse_transform(Y_pred)
    
    Weights = model.get_weights()

    result_path = './{}_DNN_{}{}{}_Results_{}.pickle'.format(scheme, hl1, hl2, hl3, j)
    with open(result_path, 'wb') as f:
        pickle.dump([Original_Testing_X, Original_Testing_Y, Prediction_Y, Weights], f)

    dest_name = './{}_DNN_{}{}{}_Results_{}.mat'.format(scheme, hl1, hl2, hl3, j)
    a = pickle.load(open(result_path, "rb"))
    scipy.io.savemat(dest_name, {
        '{}_DNN_{}{}{}_test_x_{}'.format(scheme, hl1, hl2, hl3, j): a[0],
        '{}_DNN_{}{}{}_test_y_{}'.format(scheme, hl1, hl2, hl3, j): a[1],
        '{}_DNN_{}{}{}_corrected_y_{}'.format(scheme, hl1, hl2, hl3, j): a[2],
        '{}_DNN_{}{}{}_weights_{}'.format(scheme, hl1, hl2, hl3, j): a[3]
    })
    print("Data successfully converted to .mat file ")
